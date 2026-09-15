// text_buffer.hpp -- the text model behind the terminal editor.
//
// A conf-file produced by Rockable can easily reach a hundred thousand lines,
// so the text is held as a vector of lines and the editor only ever renders the
// handful that are visible. Editing goes through two primitives, rawInsert()
// and rawErase(); everything else (typing, backspace, snippets, undo) is
// expressed with them, which keeps the undo log to a single record type.

#ifndef TEXT_BUFFER_HPP
#define TEXT_BUFFER_HPP

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class TextBuffer {
 public:
  struct Pos {
    int line = 0;
    int col = 0;
  };

  TextBuffer() { m_lines.push_back(std::string()); }

  // --- loading and saving --------------------------------------------------

  bool load(const std::string& filename) {
    std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
    if (!file) return false;

    std::ostringstream all;
    all << file.rdbuf();
    std::string content = all.str();

    m_lines.clear();
    std::string current;
    for (char c : content) {
      if (c == '\n') {
        if (!current.empty() && current.back() == '\r') current.pop_back();
        m_lines.push_back(current);
        current.clear();
      } else {
        current.push_back(c);
      }
    }
    // A trailing newline just closes the last line; anything else is a final
    // line without a newline, which we keep as such when saving.
    m_endsWithNewline = content.empty() || content.back() == '\n';
    if (!current.empty() || m_lines.empty()) m_lines.push_back(current);

    m_undo.clear();
    m_redo.clear();
    m_savedDepth = 0;
    m_cursor = Pos();
    m_selecting = false;
    return true;
  }

  bool save(const std::string& filename) {
    std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!file) return false;
    for (size_t i = 0; i < m_lines.size(); ++i) {
      file << m_lines[i];
      if (i + 1 < m_lines.size() || m_endsWithNewline) file << '\n';
    }
    if (!file) return false;
    file.close();
    m_savedDepth = m_undo.size();
    return true;
  }

  // --- inspection ----------------------------------------------------------

  int lineCount() const { return static_cast<int>(m_lines.size()); }
  const std::string& line(int i) const { return m_lines[static_cast<size_t>(i)]; }
  bool modified() const { return m_undo.size() != m_savedDepth; }
  const Pos& cursor() const { return m_cursor; }
  bool canUndo() const { return !m_undo.empty(); }
  bool canRedo() const { return !m_redo.empty(); }

  // Length in bytes; the conf files are ASCII, so this doubles as a column count.
  size_t byteSize() const {
    size_t n = 0;
    for (const std::string& l : m_lines) n += l.size() + 1;
    return n;
  }

  // The whitespace-delimited token the cursor sits on or just after, which is
  // what the documentation pane looks up.
  std::string wordAtCursor() const {
    const std::string& s = m_lines[static_cast<size_t>(m_cursor.line)];
    if (s.empty()) return std::string();
    int c = std::min(m_cursor.col, static_cast<int>(s.size()) - 1);
    if (isBlank(s[static_cast<size_t>(c)])) {
      // Standing just past the end of a token still documents that token.
      if (c == 0 || isBlank(s[static_cast<size_t>(c - 1)])) return std::string();
      c--;
    }
    int b = c;
    while (b > 0 && !isBlank(s[static_cast<size_t>(b - 1)])) b--;
    int e = c;
    while (e + 1 < static_cast<int>(s.size()) && !isBlank(s[static_cast<size_t>(e + 1)])) e++;
    return s.substr(static_cast<size_t>(b), static_cast<size_t>(e - b + 1));
  }

  // --- selection -----------------------------------------------------------
  //
  // A selection is the span between an anchor and the cursor. The editor drops
  // it on a bare movement and extends it on a movement made with Shift, so the
  // movement methods below deliberately leave it alone.

  void beginSelection() {
    if (!m_selecting) {
      m_anchor = m_cursor;
      m_selecting = true;
    }
  }

  void clearSelection() { m_selecting = false; }

  // An anchor sitting on the cursor selects nothing.
  bool hasSelection() const {
    return m_selecting && !(m_anchor.line == m_cursor.line && m_anchor.col == m_cursor.col);
  }

  Pos selectionStart() const { return isBefore(m_anchor, m_cursor) ? m_anchor : m_cursor; }
  Pos selectionEnd() const { return isBefore(m_anchor, m_cursor) ? m_cursor : m_anchor; }

  std::string selectedText() const {
    if (!hasSelection()) return std::string();
    return textInRange(selectionStart(), selectionEnd());
  }

  void deleteSelection() {
    if (!hasSelection()) return;
    const Pos from = selectionStart();
    const Pos to = selectionEnd();
    recordErase(from, to);
  }

  // --- cursor movement -----------------------------------------------------

  void setCursor(int l, int c) {
    m_cursor.line = clamp(l, 0, lineCount() - 1);
    m_cursor.col = clamp(c, 0, static_cast<int>(m_lines[static_cast<size_t>(m_cursor.line)].size()));
    m_wantCol = m_cursor.col;
  }

  void moveLeft() {
    if (m_cursor.col > 0) {
      m_cursor.col--;
    } else if (m_cursor.line > 0) {
      m_cursor.line--;
      m_cursor.col = static_cast<int>(m_lines[static_cast<size_t>(m_cursor.line)].size());
    }
    m_wantCol = m_cursor.col;
  }

  void moveRight() {
    if (m_cursor.col < static_cast<int>(m_lines[static_cast<size_t>(m_cursor.line)].size())) {
      m_cursor.col++;
    } else if (m_cursor.line + 1 < lineCount()) {
      m_cursor.line++;
      m_cursor.col = 0;
    }
    m_wantCol = m_cursor.col;
  }

  // Vertical moves remember the column the user aimed for, so that crossing a
  // short line does not lose the horizontal position.
  void moveVertical(int delta) {
    int target = clamp(m_cursor.line + delta, 0, lineCount() - 1);
    m_cursor.line = target;
    m_cursor.col = std::min(m_wantCol, static_cast<int>(m_lines[static_cast<size_t>(target)].size()));
  }

  void moveHome() {
    // First press goes to the first non-blank character, second to column 0.
    const std::string& s = m_lines[static_cast<size_t>(m_cursor.line)];
    int indent = 0;
    while (indent < static_cast<int>(s.size()) && isBlank(s[static_cast<size_t>(indent)])) indent++;
    m_cursor.col = (m_cursor.col == indent) ? 0 : indent;
    m_wantCol = m_cursor.col;
  }

  void moveEnd() {
    m_cursor.col = static_cast<int>(m_lines[static_cast<size_t>(m_cursor.line)].size());
    m_wantCol = m_cursor.col;
  }

  void moveDocStart() { setCursor(0, 0); }
  void moveDocEnd() { setCursor(lineCount() - 1, static_cast<int>(m_lines.back().size())); }

  // --- editing -------------------------------------------------------------

  void insertText(const std::string& text) {
    if (text.empty()) return;
    clearRedo();
    clearSelection();
    Pos start = m_cursor;
    Pos end = rawInsert(start, text);

    bool coalesced = false;
    if (text.find('\n') == std::string::npos && !m_undo.empty()) {
      Edit& last = m_undo.back();
      if (last.isInsert && last.coalescable && last.pos.line == start.line &&
          last.pos.col + static_cast<int>(last.text.size()) == start.col) {
        last.text += text;
        coalesced = true;
      }
    }
    if (!coalesced) {
      Edit e;
      e.isInsert = true;
      e.pos = start;
      e.text = text;
      e.cursorBefore = start;
      e.coalescable = (text.size() == 1 && text[0] != '\n');
      m_undo.push_back(e);
    }
    m_cursor = end;
    m_wantCol = end.col;
  }

  void backspace() {
    if (m_cursor.line == 0 && m_cursor.col == 0) return;
    clearRedo();
    clearSelection();
    Pos before = m_cursor;
    Pos from = before;
    if (from.col > 0) {
      from.col--;
    } else {
      from.line--;
      from.col = static_cast<int>(m_lines[static_cast<size_t>(from.line)].size());
    }
    std::string removed = rawErase(from, before);

    bool coalesced = false;
    if (removed != "\n" && !m_undo.empty()) {
      Edit& last = m_undo.back();
      if (!last.isInsert && last.coalescable && last.pos.line == from.line && last.pos.col == before.col) {
        last.text = removed + last.text;
        last.pos = from;
        coalesced = true;
      }
    }
    if (!coalesced) {
      Edit e;
      e.isInsert = false;
      e.pos = from;
      e.text = removed;
      e.cursorBefore = before;
      e.coalescable = (removed.size() == 1 && removed[0] != '\n');
      m_undo.push_back(e);
    }
    m_cursor = from;
    m_wantCol = from.col;
  }

  void deleteForward() {
    Pos from = m_cursor;
    Pos to = from;
    if (to.col < static_cast<int>(m_lines[static_cast<size_t>(to.line)].size())) {
      to.col++;
    } else if (to.line + 1 < lineCount()) {
      to.line++;
      to.col = 0;
    } else {
      return;
    }
    recordErase(from, to);
  }

  // Removes the current line and returns it (with its newline), nano style.
  std::string cutLine() {
    Pos from{m_cursor.line, 0};
    Pos to;
    if (m_cursor.line + 1 < lineCount()) {
      to = Pos{m_cursor.line + 1, 0};
    } else {
      to = Pos{m_cursor.line, static_cast<int>(m_lines[static_cast<size_t>(m_cursor.line)].size())};
    }
    if (from.line == to.line && from.col == to.col) return std::string();
    std::string text = m_lines[static_cast<size_t>(m_cursor.line)] + "\n";
    recordErase(from, to);
    return text;
  }

  void undo() { step(m_undo, m_redo); }
  void redo() { step(m_redo, m_undo); }

  // --- search --------------------------------------------------------------

  // Searches forward from just after the cursor, wrapping around once. Moves
  // the cursor onto the match and returns true when one is found.
  bool findNext(const std::string& needle, bool caseSensitive) {
    if (needle.empty()) return false;
    clearSelection();
    int startLine = m_cursor.line;
    int startCol = m_cursor.col + 1;
    for (int n = 0; n <= lineCount(); ++n) {
      int l = (startLine + n) % lineCount();
      int from = (n == 0) ? std::min(startCol, static_cast<int>(m_lines[static_cast<size_t>(l)].size())) : 0;
      size_t hit = search(m_lines[static_cast<size_t>(l)], needle, static_cast<size_t>(from), caseSensitive);
      if (hit != std::string::npos) {
        setCursor(l, static_cast<int>(hit));
        return true;
      }
    }
    return false;
  }

 private:
  struct Edit {
    bool isInsert = true;
    Pos pos;             // where the text was inserted, or removed from
    std::string text;    // the text inserted or removed
    Pos cursorBefore;    // so that undo puts the caret back where it was
    bool coalescable = false;
  };

  // Dropping a redo branch that contained the last saved state makes that
  // state unreachable, so the buffer must stay flagged as modified.
  void clearRedo() {
    if (!m_redo.empty() && m_savedDepth > m_undo.size()) m_savedDepth = static_cast<size_t>(-1);
    m_redo.clear();
  }

  static bool isBlank(char c) { return c == ' ' || c == '\t'; }

  static bool isBefore(Pos a, Pos b) {
    return a.line < b.line || (a.line == b.line && a.col < b.col);
  }

  std::string textInRange(Pos from, Pos to) const {
    if (from.line == to.line) {
      return m_lines[static_cast<size_t>(from.line)].substr(static_cast<size_t>(from.col),
                                                            static_cast<size_t>(to.col - from.col));
    }
    std::string out = m_lines[static_cast<size_t>(from.line)].substr(static_cast<size_t>(from.col));
    out += '\n';
    for (int l = from.line + 1; l < to.line; ++l) {
      out += m_lines[static_cast<size_t>(l)];
      out += '\n';
    }
    out += m_lines[static_cast<size_t>(to.line)].substr(0, static_cast<size_t>(to.col));
    return out;
  }
  static int clamp(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }

  static size_t search(const std::string& hay, const std::string& needle, size_t from, bool caseSensitive) {
    if (caseSensitive) return hay.find(needle, from);
    std::string h = hay, n = needle;
    for (char& c : h) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    for (char& c : n) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return h.find(n, from);
  }

  static std::vector<std::string> splitLines(const std::string& text) {
    std::vector<std::string> parts(1);
    for (char c : text) {
      if (c == '\n') {
        parts.push_back(std::string());
      } else {
        parts.back().push_back(c);
      }
    }
    return parts;
  }

  Pos rawInsert(Pos at, const std::string& text) {
    std::vector<std::string> parts = splitLines(text);
    std::string& target = m_lines[static_cast<size_t>(at.line)];
    std::string tail = target.substr(static_cast<size_t>(at.col));
    target.erase(static_cast<size_t>(at.col));
    target += parts[0];
    if (parts.size() == 1) {
      target += tail;
      return Pos{at.line, at.col + static_cast<int>(parts[0].size())};
    }
    std::vector<std::string> extra(parts.begin() + 1, parts.end());
    Pos end{at.line + static_cast<int>(extra.size()), static_cast<int>(extra.back().size())};
    extra.back() += tail;
    m_lines.insert(m_lines.begin() + at.line + 1, extra.begin(), extra.end());
    return end;
  }

  std::string rawErase(Pos from, Pos to) {
    const std::string removed = textInRange(from, to);
    if (from.line == to.line) {
      m_lines[static_cast<size_t>(from.line)].erase(static_cast<size_t>(from.col),
                                                    static_cast<size_t>(to.col - from.col));
      return removed;
    }
    std::string merged = m_lines[static_cast<size_t>(from.line)].substr(0, static_cast<size_t>(from.col)) +
                         m_lines[static_cast<size_t>(to.line)].substr(static_cast<size_t>(to.col));
    m_lines.erase(m_lines.begin() + from.line, m_lines.begin() + to.line + 1);
    m_lines.insert(m_lines.begin() + from.line, merged);
    return removed;
  }

  void recordErase(Pos from, Pos to) {
    clearRedo();
    clearSelection();
    Edit e;
    e.isInsert = false;
    e.pos = from;
    e.cursorBefore = m_cursor;
    e.text = rawErase(from, to);
    e.coalescable = false;
    m_undo.push_back(e);
    m_cursor = from;
    m_wantCol = from.col;
  }

  static Pos advance(Pos from, const std::string& text) {
    Pos p = from;
    for (char c : text) {
      if (c == '\n') {
        p.line++;
        p.col = 0;
      } else {
        p.col++;
      }
    }
    return p;
  }

  // Undo and redo are the same operation seen from either side: pop the last
  // edit off one stack, apply its inverse, and push the inverse on the other.
  void step(std::vector<Edit>& from, std::vector<Edit>& to) {
    if (from.empty()) return;
    clearSelection();
    Edit e = from.back();
    from.pop_back();

    Edit inverse = e;
    inverse.isInsert = !e.isInsert;
    inverse.coalescable = false;
    inverse.cursorBefore = m_cursor;  // where to put the caret back when replayed

    if (e.isInsert) {
      rawErase(e.pos, advance(e.pos, e.text));
    } else {
      rawInsert(e.pos, e.text);
    }
    m_cursor = e.cursorBefore;
    m_cursor.line = clamp(m_cursor.line, 0, lineCount() - 1);
    m_cursor.col = clamp(m_cursor.col, 0, static_cast<int>(m_lines[static_cast<size_t>(m_cursor.line)].size()));
    m_wantCol = m_cursor.col;
    to.push_back(inverse);
  }

  std::vector<std::string> m_lines;
  Pos m_cursor;
  Pos m_anchor;             // the other end of the selection
  bool m_selecting = false;
  int m_wantCol = 0;
  size_t m_savedDepth = 0;  // undo depth at the last save; drives modified()
  bool m_endsWithNewline = true;
  std::vector<Edit> m_undo;
  std::vector<Edit> m_redo;
};

#endif /* end of include guard: TEXT_BUFFER_HPP */
