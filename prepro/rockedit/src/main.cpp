// rockedit -- a small terminal editor for the Rockable input files.
//
// It is the terminal counterpart of prepro/confedit: same syntax colouring and
// same inline documentation, but light enough to be used over ssh on a cluster.
// Both editors read their keywords, types, documentation and snippets from
// prepro/common/rockable.lang, so the two cannot drift apart.
//
// Usage: rockedit [file]
// Press ^G once inside for the key map.

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "rockable_lang.hpp"
#include "terminal.hpp"
#include "text_buffer.hpp"

namespace {

// The four styles of confedit: plain text, comment, keyword, enumerated value.
enum class Style : uint8_t { Plain, Comment, Keyword, Type };

std::string styleCode(Style style) {
  switch (style) {
    case Style::Comment: return term::sgr::green;
    case Style::Keyword: return term::sgr::blue;
    case Style::Type: return term::sgr::cyan;
    default: return std::string();
  }
}

// Same rule as confedit's style_parse(): split the line on whitespace, treat a
// token starting with '/', '#' or '!' as a comment running to the end of the
// line, and colour the remaining tokens by looking them up in the language.
std::vector<Style> styleLine(const std::string& s, const rockable_lang::Language& lang) {
  std::vector<Style> style(s.size(), Style::Plain);
  size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) i++;
    if (i >= s.size()) break;
    const size_t begin = i;
    while (i < s.size() && !std::isspace(static_cast<unsigned char>(s[i]))) i++;

    if (s[begin] == '/' || s[begin] == '#' || s[begin] == '!') {
      for (size_t k = begin; k < s.size(); ++k) style[k] = Style::Comment;
      return style;
    }
    switch (lang.kindOf(s.substr(begin, i - begin))) {
      case rockable_lang::Kind::Keyword:
        for (size_t k = begin; k < i; ++k) style[k] = Style::Keyword;
        break;
      case rockable_lang::Kind::Type:
        for (size_t k = begin; k < i; ++k) style[k] = Style::Type;
        break;
      default:
        break;
    }
  }
  return style;
}

// Greedy word wrap, so that the documentation pane can be given the exact
// height it needs.
std::vector<std::string> wrap(const std::string& text, int width) {
  std::vector<std::string> out;
  if (width < 8) width = 8;
  size_t begin = 0;
  while (begin <= text.size()) {
    const size_t eol = text.find('\n', begin);
    std::string paragraph = text.substr(begin, (eol == std::string::npos) ? std::string::npos : eol - begin);

    std::string indent;  // continuation lines keep the paragraph's indentation
    for (char c : paragraph) {
      if (c != ' ') break;
      indent += ' ';
    }
    if (static_cast<int>(indent.size()) > width / 2) indent.clear();

    while (static_cast<int>(paragraph.size()) > width) {
      size_t cut = paragraph.rfind(' ', static_cast<size_t>(width));
      if (cut == std::string::npos || cut <= indent.size()) cut = static_cast<size_t>(width);
      out.push_back(paragraph.substr(0, cut));
      const size_t next = paragraph.find_first_not_of(' ', cut);
      if (next == std::string::npos) {
        paragraph.clear();
        break;
      }
      paragraph = indent + paragraph.substr(next);
    }
    out.push_back(paragraph);
    if (eol == std::string::npos) break;
    begin = eol + 1;
  }
  return out;
}

std::string basename(const std::string& path) {
  const size_t slash = path.find_last_of('/');
  return (slash == std::string::npos) ? path : path.substr(slash + 1);
}

std::string lowered(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

class Editor {
 public:
  Editor(term::Terminal* terminal, rockable_lang::Language* lang, std::string filename)
      : m_terminal(terminal), m_lang(lang), m_filename(std::move(filename)) {}

  TextBuffer& buffer() { return m_buffer; }
  void setMessage(std::string message) { m_message = std::move(message); }

  void run() {
    while (!m_quit) {
      draw();
      const int key = m_terminal->readKey();
      if (key == term::KeyNone) continue;
      if (key == term::KeyResize) continue;  // the next draw() picks up the new size
      onKey(key);
    }
  }

 private:
  enum class Mode { Normal, Prompt, Snippets, Help };
  enum class Prompt { Find, GotoLine, SaveAs, QuitConfirm };

  // --- drawing -------------------------------------------------------------

  void draw() {
    const term::Size size = m_terminal->size();
    const int width = std::max(24, size.cols);
    const int gutter = gutterWidth();
    const int contentWidth = std::max(8, width - gutter);

    std::vector<std::string> docRows;
    if (m_showDoc) docRows = documentationRows(width);
    // header + status + shortcuts, plus the documentation pane and its rule
    const int chrome = 3 + (m_showDoc ? static_cast<int>(docRows.size()) + 1 : 0);
    const int textHeight = std::max(1, size.rows - chrome);

    scrollIntoView(textHeight, contentWidth);

    term::Frame frame(term::Size{size.rows, width});
    frame.add(headerRow(width));

    for (int i = 0; i < textHeight; ++i) {
      const int l = m_top + i;
      if (l >= m_buffer.lineCount()) {
        term::Row row(width);
        row.put("~", term::sgr::dim);
        frame.add(row);
      } else {
        frame.add(textRow(l, gutter, contentWidth, width));
      }
    }

    if (m_showDoc) {
      term::Row rule(width);
      rule.put(std::string(static_cast<size_t>(width), '-'), term::sgr::dim);
      frame.add(rule);
      for (const std::string& line : docRows) {
        term::Row row(width);
        row.put(line);
        frame.add(row);
      }
    }

    frame.add(statusRow(width));
    frame.add(shortcutRow(width));

    // The caret sits either in the text or at the end of the prompt.
    int caretRow = 1 + (m_buffer.cursor().line - m_top);
    int caretCol = gutter + (m_buffer.cursor().col - m_left);
    if (m_mode == Mode::Prompt) {
      caretRow = size.rows - 2;
      caretCol = std::min(width - 1, static_cast<int>(m_promptLabel.size() + m_promptText.size()));
    }

    if (m_mode == Mode::Snippets) {
      drawSnippetOverlay(frame, size, width);
      return;
    }
    if (m_mode == Mode::Help) {
      drawHelpOverlay(frame, size, width);
      return;
    }
    frame.flush(caretRow, caretCol, true);
  }

  int gutterWidth() const {
    int digits = 1;
    for (int n = m_buffer.lineCount(); n >= 10; n /= 10) digits++;
    return digits + 2;  // the number, a space, and the gutter separator
  }

  void scrollIntoView(int textHeight, int contentWidth) {
    const TextBuffer::Pos c = m_buffer.cursor();
    if (c.line < m_top) m_top = c.line;
    if (c.line >= m_top + textHeight) m_top = c.line - textHeight + 1;
    m_top = std::max(0, std::min(m_top, std::max(0, m_buffer.lineCount() - 1)));

    if (c.col < m_left) m_left = c.col;
    if (c.col >= m_left + contentWidth) m_left = c.col - contentWidth + 1;
    m_left = std::max(0, m_left);
  }

  term::Row textRow(int l, int gutter, int contentWidth, int width) const {
    const std::string& s = m_buffer.line(l);
    const std::vector<Style> style = styleLine(s, *m_lang);
    const bool isCursorLine = (l == m_buffer.cursor().line);

    // The part of this line covered by the selection, if any. A selected line
    // break shows as one highlighted column past the end of the text.
    int selectionFrom = -1;
    int selectionTo = -1;
    if (m_buffer.hasSelection()) {
      const TextBuffer::Pos begin = m_buffer.selectionStart();
      const TextBuffer::Pos end = m_buffer.selectionEnd();
      if (l >= begin.line && l <= end.line) {
        selectionFrom = (l == begin.line) ? begin.col : 0;
        selectionTo = (l == end.line) ? end.col : static_cast<int>(s.size()) + 1;
      }
    }

    term::Row row(width);
    std::string number = std::to_string(l + 1);
    if (static_cast<int>(number.size()) < gutter - 1) {
      number.insert(0, static_cast<size_t>(gutter - 1) - number.size(), ' ');
    }
    number += ' ';
    row.put(number, isCursorLine ? term::sgr::yellow : term::sgr::dim);

    // Emit runs of identical style rather than one escape sequence per byte.
    const int last = std::min(m_left + contentWidth, std::max(static_cast<int>(s.size()), selectionTo));
    std::string run;
    Style runStyle = Style::Plain;
    bool runSelected = false;
    auto flush = [&]() {
      if (run.empty()) return;
      row.put(run, (runSelected ? std::string(term::sgr::inverse) : std::string()) + styleCode(runStyle));
      run.clear();
    };

    for (int c = m_left; c < last; ++c) {
      const bool inLine = c < static_cast<int>(s.size());
      char ch = inLine ? s[static_cast<size_t>(c)] : ' ';
      if (ch == '\t') ch = ' ';  // tabs are unusual in conf files; show one space
      const Style st = inLine ? style[static_cast<size_t>(c)] : Style::Plain;
      const bool selected = (selectionFrom >= 0 && c >= selectionFrom && c < selectionTo);

      if (!run.empty() && (st != runStyle || selected != runSelected)) flush();
      if (run.empty()) {
        runStyle = st;
        runSelected = selected;
      }
      run.push_back(ch);
    }
    flush();
    return row;
  }

  term::Row headerRow(int width) const {
    term::Row row(width);
    row.put(" rockedit ", term::sgr::inverse);
    row.put("  " + (m_filename.empty() ? std::string("[new file]") : basename(m_filename)));
    if (m_buffer.modified()) row.put("  *modified*", term::sgr::yellow);

    const std::string right = m_lang->loaded()
                                  ? (std::to_string(m_lang->wordCount()) + " keywords ")
                                  : std::string("rockable.lang NOT found ");
    row.padTo(width - static_cast<int>(right.size()));
    row.put(right, m_lang->loaded() ? term::sgr::dim : term::sgr::red);
    return row;
  }

  term::Row statusRow(int width) const {
    term::Row row(width);
    if (m_mode == Mode::Prompt) {
      row.put(m_promptLabel, term::sgr::yellow);
      row.put(m_promptText);
      return row;
    }
    const TextBuffer::Pos c = m_buffer.cursor();
    row.put(" Ln " + std::to_string(c.line + 1) + "/" + std::to_string(m_buffer.lineCount()) + ", Col " +
                std::to_string(c.col + 1),
            term::sgr::dim);
    if (m_buffer.hasSelection()) {
      const TextBuffer::Pos begin = m_buffer.selectionStart();
      const TextBuffer::Pos end = m_buffer.selectionEnd();
      const int selectedLines = end.line - begin.line + 1;
      row.put(selectedLines > 1 ? ("   [" + std::to_string(selectedLines) + " lines selected]")
                                : ("   [" + std::to_string(end.col - begin.col) + " chars selected]"),
              term::sgr::yellow);
    }
    if (!m_message.empty()) {
      row.put("   " + m_message, term::sgr::yellow);
    } else if (!m_showDoc && !m_buffer.hasSelection()) {
      const std::string word = documentedWord();
      if (!word.empty() && m_lang->doc(word) != nullptr) {
        row.put("   ^D documents '" + word + "'", term::sgr::dim);
      }
    }
    return row;
  }

  term::Row shortcutRow(int width) const {
    struct Key {
      const char* key;
      const char* name;
    };
    static const Key keys[] = {{"^O", "Save"}, {"^W", "Find"},  {"^N", "Next"},  {"^D", "Doc"},
                               {"^P", "Snip"}, {"^L", "GoTo"},  {"^K", "Cut"},   {"^U", "Paste"},
                               {"^Z", "Undo"}, {"^G", "Help"},  {"^X", "Quit"}};
    term::Row row(width);
    for (const Key& k : keys) {
      if (row.room() < 8) break;
      row.put(k.key, term::sgr::inverse);
      row.put(std::string(" ") + k.name + " ");
    }
    return row;
  }

  // confedit documents the selected text; with no selection we fall back to the
  // token under the cursor.
  std::string documentedWord() const {
    if (!m_buffer.hasSelection()) return m_buffer.wordAtCursor();
    std::string selected = m_buffer.selectedText();
    const size_t begin = selected.find_first_not_of(" \t\n");
    if (begin == std::string::npos) return std::string();
    const size_t end = selected.find_last_not_of(" \t\n");
    return selected.substr(begin, end - begin + 1);
  }

  std::vector<std::string> documentationRows(int width) const {
    const std::string word = documentedWord();
    if (word.empty()) return {" no word under the cursor"};

    const std::string* doc = m_lang->doc(word);
    if (doc == nullptr) return {" " + word + "  is not documented"};

    std::vector<std::string> rows;
    rows.push_back(" " + word);
    const std::vector<std::string> lines = wrap(*doc, width - 2);
    const size_t maxRows = 12;
    for (size_t i = 0; i < lines.size() && i < maxRows; ++i) rows.push_back(" " + lines[i]);
    if (lines.size() > maxRows) rows.push_back(" ...");
    return rows;
  }

  // The overlays repaint the middle of the screen on top of a finished frame,
  // which keeps them independent of the editor layout.
  void drawSnippetOverlay(term::Frame& frame, term::Size size, int width) {
    frame.flush(0, 0, false);

    const std::vector<size_t> matches = matchingSnippets();
    const int boxWidth = std::min(std::max(44, width - 8), 78);
    const int boxLeft = std::max(0, (width - boxWidth) / 2);
    const int listRows = std::max(3, std::min(10, size.rows - 12));

    std::vector<std::string> lines;
    lines.push_back(" Snippets");
    lines.push_back(" filter: " + m_snippetFilter + "_");
    lines.push_back(std::string(static_cast<size_t>(boxWidth - 2), '-'));

    int first = 0;
    if (m_snippetSel >= listRows) first = m_snippetSel - listRows + 1;
    if (matches.empty()) {
      lines.push_back(" (no match)");
    } else {
      for (int i = first; i < first + listRows && i < static_cast<int>(matches.size()); ++i) {
        const std::string& title = m_lang->snippets()[matches[static_cast<size_t>(i)]].title;
        lines.push_back((i == m_snippetSel ? " > " : "   ") + title);
      }
    }
    if (!matches.empty() && m_snippetSel < static_cast<int>(matches.size())) {
      lines.push_back(std::string(static_cast<size_t>(boxWidth - 2), '-'));
      const std::string& body = m_lang->snippets()[matches[static_cast<size_t>(m_snippetSel)]].body;
      for (const std::string& l : wrap(body, boxWidth - 4)) lines.push_back("  " + l);
    }
    lines.push_back(std::string(static_cast<size_t>(boxWidth - 2), '-'));
    lines.push_back(" Enter insert   Esc cancel");

    paintBox(lines, boxLeft, boxWidth, size);
  }

  void drawHelpOverlay(term::Frame& frame, term::Size size, int width) {
    frame.flush(0, 0, false);

    struct Row {
      const char* key;
      const char* what;
    };
    static const Row help[] = {
        {"^O or ^S", "save (asks for a name when there is none)"},
        {"^X", "quit"},
        {"^W", "find; ^N finds the next occurrence (case-insensitive)"},
        {"^L", "go to a line number"},
        {"^D", "documentation of the selection, or of the word under the cursor"},
        {"^P", "insert a snippet"},
        {"Shift+arrows", "select; Shift+Home/End/PgUp/PgDn too"},
        {"^K", "cut the selection, or the line when there is none; ^U pastes"},
        {"^Z", "undo; ^Y redo"},
        {"^A / ^E", "start / end of the line (Home and End too)"},
        {"Ctrl+Up/Dn", "start / end of the file"},
        {"^G", "this help"},
    };

    std::vector<std::string> lines;
    lines.push_back(" rockedit -- keys");
    lines.push_back("");
    for (const Row& r : help) {
      std::string key = r.key;
      key.resize(14, ' ');
      lines.push_back("  " + key + r.what);
    }
    lines.push_back("");
    lines.push_back(m_lang->loaded() ? ("  language: " + m_lang->path())
                                     : std::string("  rockable.lang was NOT found"));
    lines.push_back("");
    lines.push_back("  Press any key to close.");

    size_t boxWidth = 40;
    for (const std::string& l : lines) boxWidth = std::max(boxWidth, l.size() + 3);
    boxWidth = std::min(boxWidth, static_cast<size_t>(std::max(40, width - 4)));
    paintBox(lines, std::max(0, (width - static_cast<int>(boxWidth)) / 2), static_cast<int>(boxWidth), size);
  }

  // Draws a framed box over the current screen contents.
  static void paintBox(const std::vector<std::string>& lines, int left, int boxWidth, term::Size size) {
    const int boxHeight = static_cast<int>(lines.size()) + 2;
    const int top = std::max(0, (size.rows - boxHeight) / 2);
    std::string out;
    char move[32];

    auto emit = [&](int row, const std::string& content) {
      if (row < 0 || row >= size.rows) return;
      snprintf(move, sizeof(move), "\x1b[%d;%dH", row + 1, left + 1);
      out += move;
      out += content;
    };

    const std::string horizontal(static_cast<size_t>(boxWidth - 2), '-');
    emit(top, std::string(term::sgr::inverse) + "+" + horizontal + "+" + term::sgr::reset);
    for (size_t i = 0; i < lines.size(); ++i) {
      std::string content = lines[i];
      if (static_cast<int>(content.size()) > boxWidth - 2) content.resize(static_cast<size_t>(boxWidth - 2));
      content.resize(static_cast<size_t>(boxWidth - 2), ' ');
      emit(top + 1 + static_cast<int>(i),
           std::string(term::sgr::inverse) + "|" + term::sgr::reset + content + term::sgr::inverse + "|" +
               term::sgr::reset);
    }
    emit(top + boxHeight - 1, std::string(term::sgr::inverse) + "+" + horizontal + "+" + term::sgr::reset);
    out += "\x1b[?25l";
    term::Terminal::put(out);
  }

  std::vector<size_t> matchingSnippets() const {
    std::vector<size_t> out;
    const std::string needle = lowered(m_snippetFilter);
    for (size_t i = 0; i < m_lang->snippets().size(); ++i) {
      if (needle.empty() || lowered(m_lang->snippets()[i].title).find(needle) != std::string::npos) {
        out.push_back(i);
      }
    }
    return out;
  }

  // --- keys ----------------------------------------------------------------

  void onKey(int key) {
    m_message.clear();
    switch (m_mode) {
      case Mode::Help: m_mode = Mode::Normal; return;
      case Mode::Snippets: onSnippetKey(key); return;
      case Mode::Prompt: onPromptKey(key); return;
      case Mode::Normal: onNormalKey(key); return;
    }
  }

  static bool isPrintable(int key) { return key >= 32 && key < 127; }

  void onSnippetKey(int key) {
    const std::vector<size_t> matches = matchingSnippets();
    if (key == term::KeyEscape) {
      m_mode = Mode::Normal;
      return;
    }
    if (key == term::KeyArrowDown) {
      if (m_snippetSel + 1 < static_cast<int>(matches.size())) m_snippetSel++;
      return;
    }
    if (key == term::KeyArrowUp) {
      if (m_snippetSel > 0) m_snippetSel--;
      return;
    }
    if (key == term::KeyReturn) {
      if (!matches.empty() && m_snippetSel < static_cast<int>(matches.size())) {
        const rockable_lang::Snippet& s = m_lang->snippets()[matches[static_cast<size_t>(m_snippetSel)]];
        m_buffer.moveEnd();
        m_buffer.insertText("\n" + s.body);
        m_message = "inserted snippet '" + s.title + "'";
      }
      m_mode = Mode::Normal;
      return;
    }
    if (key == term::KeyBackspace || key == term::ctrl('H')) {
      if (!m_snippetFilter.empty()) m_snippetFilter.pop_back();
      m_snippetSel = 0;
      return;
    }
    if (isPrintable(key)) {
      m_snippetFilter.push_back(static_cast<char>(key));
      m_snippetSel = 0;
    }
  }

  void onPromptKey(int key) {
    if (key == term::KeyEscape) {
      m_mode = Mode::Normal;
      m_message = "cancelled";
      return;
    }

    if (m_prompt == Prompt::QuitConfirm) {
      if (key == 'y' || key == 'Y') {
        m_mode = Mode::Normal;
        if (saveToCurrentName()) m_quit = true;
        return;
      }
      if (key == 'n' || key == 'N') m_quit = true;
      return;
    }

    if (key == term::KeyBackspace || key == term::ctrl('H')) {
      if (!m_promptText.empty()) m_promptText.pop_back();
      return;
    }
    if (key == term::KeyReturn) {
      const std::string answer = m_promptText;
      const Prompt kind = m_prompt;
      m_mode = Mode::Normal;
      applyPrompt(kind, answer);
      return;
    }
    if (isPrintable(key)) m_promptText.push_back(static_cast<char>(key));
  }

  void applyPrompt(Prompt kind, const std::string& answer) {
    switch (kind) {
      case Prompt::Find:
        if (answer.empty()) return;
        m_lastSearch = answer;
        if (!m_buffer.findNext(m_lastSearch, false)) m_message = "'" + answer + "' not found";
        return;

      case Prompt::GotoLine: {
        const long n = std::strtol(answer.c_str(), nullptr, 10);
        if (n <= 0) {
          m_message = "not a line number";
          return;
        }
        m_buffer.setCursor(static_cast<int>(n) - 1, 0);
        return;
      }

      case Prompt::SaveAs:
        if (answer.empty()) {
          m_message = "save cancelled";
          return;
        }
        m_filename = answer;
        saveToCurrentName();
        return;

      case Prompt::QuitConfirm:
        return;
    }
  }

  bool saveToCurrentName() {
    if (m_filename.empty()) {
      openPrompt(Prompt::SaveAs, "Save as: ", "");
      return false;
    }
    if (m_buffer.save(m_filename)) {
      m_message = "written to " + m_filename;
      return true;
    }
    m_message = "CANNOT write to " + m_filename;
    return false;
  }

  void openPrompt(Prompt kind, const std::string& label, const std::string& initial) {
    m_mode = Mode::Prompt;
    m_prompt = kind;
    m_promptLabel = label;
    m_promptText = initial;
  }

  // Shift extends the selection from wherever it was anchored; a bare movement
  // drops it.
  void startMove(bool extend) {
    if (extend) {
      m_buffer.beginSelection();
    } else {
      m_buffer.clearSelection();
    }
  }

  // Typing, pasting or deleting replaces whatever is selected.
  void dropSelection() {
    if (m_buffer.hasSelection()) m_buffer.deleteSelection();
  }

  void onNormalKey(int key) {
    switch (key) {
      // --- movement, with or without Shift
      case term::KeyArrowLeft:
      case term::KeyShiftArrowLeft:
        startMove(key == term::KeyShiftArrowLeft);
        m_buffer.moveLeft();
        return;
      case term::KeyArrowRight:
      case term::KeyShiftArrowRight:
        startMove(key == term::KeyShiftArrowRight);
        m_buffer.moveRight();
        return;
      case term::KeyArrowUp:
      case term::KeyShiftArrowUp:
        startMove(key == term::KeyShiftArrowUp);
        m_buffer.moveVertical(-1);
        return;
      case term::KeyArrowDown:
      case term::KeyShiftArrowDown:
        startMove(key == term::KeyShiftArrowDown);
        m_buffer.moveVertical(1);
        return;
      case term::KeyPageUp:
      case term::KeyShiftPageUp:
        startMove(key == term::KeyShiftPageUp);
        m_buffer.moveVertical(-pageSize());
        return;
      case term::KeyPageDown:
      case term::KeyShiftPageDown:
        startMove(key == term::KeyShiftPageDown);
        m_buffer.moveVertical(pageSize());
        return;
      case term::KeyHome:
      case term::KeyShiftHome:
        startMove(key == term::KeyShiftHome);
        m_buffer.moveHome();
        return;
      case term::KeyEnd:
      case term::KeyShiftEnd:
        startMove(key == term::KeyShiftEnd);
        m_buffer.moveEnd();
        return;
      case term::KeyCtrlArrowUp:
        startMove(false);
        m_buffer.moveDocStart();
        return;
      case term::KeyCtrlArrowDown:
        startMove(false);
        m_buffer.moveDocEnd();
        return;

      // --- editing
      case term::KeyReturn:
        dropSelection();
        m_buffer.insertText("\n");
        return;
      case term::KeyTab:
        dropSelection();
        m_buffer.insertText("  ");
        return;
      case term::KeyBackspace:
        if (m_buffer.hasSelection()) {
          m_buffer.deleteSelection();
        } else {
          m_buffer.backspace();
        }
        return;
      case term::KeyDelete:
        if (m_buffer.hasSelection()) {
          m_buffer.deleteSelection();
        } else {
          m_buffer.deleteForward();
        }
        return;

      default: break;
    }

    if (key == term::ctrl('A')) { startMove(false); m_buffer.moveHome(); return; }
    if (key == term::ctrl('E')) { startMove(false); m_buffer.moveEnd(); return; }
    if (key == term::ctrl('H')) {
      if (m_buffer.hasSelection()) {
        m_buffer.deleteSelection();
      } else {
        m_buffer.backspace();
      }
      return;
    }

    if (key == term::ctrl('X')) {
      if (m_buffer.modified()) {
        openPrompt(Prompt::QuitConfirm, "Save the modified buffer? (y/n, Esc goes back) ", "");
      } else {
        m_quit = true;
      }
      return;
    }
    // XON/XOFF is off, so ^S is free and works as an alias for ^O.
    if (key == term::ctrl('O') || key == term::ctrl('S')) { saveToCurrentName(); return; }
    if (key == term::ctrl('C')) { m_message = "use ^X to quit"; return; }
    if (key == term::ctrl('G')) { m_mode = Mode::Help; return; }
    if (key == term::ctrl('D')) { m_showDoc = !m_showDoc; return; }
    if (key == term::ctrl('W')) { openPrompt(Prompt::Find, "Find: ", m_lastSearch); return; }
    if (key == term::ctrl('L')) { openPrompt(Prompt::GotoLine, "Go to line: ", ""); return; }
    if (key == term::ctrl('N')) {
      if (m_lastSearch.empty()) {
        openPrompt(Prompt::Find, "Find: ", "");
      } else if (!m_buffer.findNext(m_lastSearch, false)) {
        m_message = "'" + m_lastSearch + "' not found";
      }
      return;
    }
    if (key == term::ctrl('P')) {
      if (m_lang->snippets().empty()) {
        m_message = "no snippet available (rockable.lang not loaded?)";
      } else {
        m_mode = Mode::Snippets;
        m_snippetFilter.clear();
        m_snippetSel = 0;
      }
      return;
    }
    if (key == term::ctrl('K')) {
      // With a selection ^K cuts exactly it; without one it cuts the line, and
      // ^U then has to put it back as a whole line rather than in the middle of
      // another one.
      if (m_buffer.hasSelection()) {
        m_clipboard = m_buffer.selectedText();
        m_clipboardIsWholeLine = false;
        m_buffer.deleteSelection();
      } else {
        const std::string cut = m_buffer.cutLine();
        if (!cut.empty()) {
          m_clipboard = cut;
          m_clipboardIsWholeLine = true;
        }
      }
      return;
    }
    if (key == term::ctrl('U')) {
      if (m_clipboard.empty()) {
        m_message = "nothing to paste";
        return;
      }
      dropSelection();
      if (m_clipboardIsWholeLine) m_buffer.setCursor(m_buffer.cursor().line, 0);
      m_buffer.insertText(m_clipboard);
      return;
    }
    if (key == term::ctrl('Z')) { m_buffer.undo(); return; }
    if (key == term::ctrl('Y')) { m_buffer.redo(); return; }

    if (isPrintable(key)) {
      dropSelection();
      m_buffer.insertText(std::string(1, static_cast<char>(key)));
    }
  }

  int pageSize() const { return std::max(1, m_terminal->size().rows - 5); }

  term::Terminal* m_terminal;
  rockable_lang::Language* m_lang;
  TextBuffer m_buffer;
  std::string m_filename;

  int m_top = 0;   // first visible line
  int m_left = 0;  // first visible column
  bool m_showDoc = false;
  bool m_quit = false;
  std::string m_message;
  std::string m_clipboard;
  bool m_clipboardIsWholeLine = false;
  std::string m_lastSearch;

  Mode m_mode = Mode::Normal;
  Prompt m_prompt = Prompt::Find;
  std::string m_promptLabel;
  std::string m_promptText;

  std::string m_snippetFilter;
  int m_snippetSel = 0;
};

}  // namespace

int main(int argc, char** argv) {
  std::string filename;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      printf("usage: %s [file]\n\n", argv[0]);
      printf("A terminal editor for the Rockable input files, with the syntax colouring\n");
      printf("and the inline documentation of confedit. Press ^G inside for the key map.\n\n");
      printf("The language is read from rockable.lang, looked up in $ROCKABLE_LANG, the\n");
      printf("current directory, next to the executable, then in ~/.rockable.\n");
      return 0;
    }
    filename = arg;
  }

  rockable_lang::Language lang;
  const std::string langFile = rockable_lang::findFile(argv[0]);
  if (!langFile.empty()) lang.load(langFile);

  term::Terminal terminal;
  if (!terminal.enter()) {
    fprintf(stderr, "rockedit: not running on a terminal\n");
    return 1;
  }

  {
    Editor editor(&terminal, &lang, filename);
    if (!filename.empty()) {
      if (editor.buffer().load(filename)) {
        editor.setMessage("read " + filename);
      } else {
        editor.setMessage(filename + " is a new file");
      }
    }
    if (langFile.empty()) {
      editor.setMessage("rockable.lang not found: no colouring, no documentation, no snippet");
    }
    editor.run();
  }

  terminal.leave();
  return 0;
}
