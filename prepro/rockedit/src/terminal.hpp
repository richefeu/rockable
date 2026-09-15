// terminal.hpp -- the thin terminal layer of rockedit: raw mode, key decoding
// and frame output, using nothing but termios and ANSI escape sequences.
//
// Nothing here knows about Rockable. The editor draws a whole frame into a
// single string and writes it in one go, which is simple and flicker-free
// enough for a text editor, even over a slow ssh link.
//
// Bytes are assumed to be columns: the Rockable input files are ASCII.

#ifndef TERMINAL_HPP
#define TERMINAL_HPP

#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/ioctl.h>
#include <sys/select.h>
#include <sys/time.h>
#include <termios.h>
#include <unistd.h>

namespace term {

// --- SGR sequences ----------------------------------------------------------
namespace sgr {
inline const char* reset = "\x1b[0m";
inline const char* bold = "\x1b[1m";
inline const char* dim = "\x1b[90m";      // bright black, readable on both themes
inline const char* inverse = "\x1b[7m";
inline const char* red = "\x1b[31m";
inline const char* green = "\x1b[32m";     // comments, as in confedit
inline const char* yellow = "\x1b[33m";
inline const char* blue = "\x1b[94m";      // keywords (bright, legible on black)
inline const char* cyan = "\x1b[36m";      // types and enumerated values
}  // namespace sgr

// --- keys -------------------------------------------------------------------
// Codes 1..26 are Ctrl+A..Ctrl+Z as the terminal delivers them, 32..126 are
// printable characters; everything else gets a symbolic value above 255.
enum Key : int {
  KeyNone = -1,
  KeyCtrl = 0,  // helper base: KeyCtrl + 'A' - 'A' + 1 == 1
  KeyTab = 9,
  KeyReturn = 13,
  KeyEscape = 27,
  KeyBackspace = 127,

  KeyArrowUp = 1000,
  KeyArrowDown,
  KeyArrowLeft,
  KeyArrowRight,
  KeyCtrlArrowUp,
  KeyCtrlArrowDown,
  KeyCtrlArrowLeft,
  KeyCtrlArrowRight,
  KeyShiftArrowUp,
  KeyShiftArrowDown,
  KeyShiftArrowLeft,
  KeyShiftArrowRight,
  KeyHome,
  KeyEnd,
  KeyShiftHome,
  KeyShiftEnd,
  KeyPageUp,
  KeyPageDown,
  KeyShiftPageUp,
  KeyShiftPageDown,
  KeyDelete,
  KeyResize
};

inline constexpr int ctrl(char c) { return c - 'A' + 1; }

struct Size {
  int rows = 24;
  int cols = 80;
};

// Set from the SIGWINCH handler; read() then fails with EINTR and readKey()
// turns that into a KeyResize.
inline volatile sig_atomic_t g_resized = 0;
inline void onWinch(int) { g_resized = 1; }

// Takes over the terminal for the lifetime of the object: raw mode plus the
// alternate screen buffer, so the user's scrollback is left untouched.
class Terminal {
 public:
  bool enter() {
    if (!isatty(STDIN_FILENO)) return false;
    if (tcgetattr(STDIN_FILENO, &m_saved) == -1) return false;

    struct termios raw = m_saved;
    // No echo, no line buffering, no signals (so that ^C and ^Z reach us as
    // ordinary keys), no XON/XOFF (so that ^S and ^Q would be usable too), no
    // CR/LF translation on input or output.
    raw.c_iflag &= ~(unsigned)(IXON | ICRNL | BRKINT | INPCK | ISTRIP);
    raw.c_oflag &= ~(unsigned)(OPOST);
    raw.c_lflag &= ~(unsigned)(ECHO | ICANON | ISIG | IEXTEN);
    raw.c_cflag |= (unsigned)CS8;
    raw.c_cc[VMIN] = 1;
    raw.c_cc[VTIME] = 0;
    if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw) == -1) return false;

    m_active = true;
    struct sigaction sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.sa_handler = onWinch;
    sigaction(SIGWINCH, &sa, nullptr);

    put("\x1b[?1049h");  // alternate screen
    put("\x1b[2J");      // and start from a clean one
    return true;
  }

  void leave() {
    if (!m_active) return;
    put("\x1b[0m");
    put("\x1b[?25h");    // make sure the cursor is visible again
    put("\x1b[?1049l");  // back to the user's screen
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &m_saved);
    m_active = false;
  }

  ~Terminal() { leave(); }

  Size size() const {
    struct winsize ws;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 0) {
      return Size{static_cast<int>(ws.ws_row), static_cast<int>(ws.ws_col)};
    }
    return Size{};
  }

  static void put(const std::string& s) {
    ssize_t written = 0;
    while (written < static_cast<ssize_t>(s.size())) {
      ssize_t n = ::write(STDOUT_FILENO, s.data() + written, s.size() - static_cast<size_t>(written));
      if (n <= 0) {
        if (n < 0 && errno == EINTR) continue;
        return;
      }
      written += n;
    }
  }

  // Blocks until a key is available, decoding the escape sequences of the
  // usual xterm-compatible terminals.
  int readKey() {
    int c = readByte(-1);
    if (c == KeyNone) return g_resized ? (g_resized = 0, KeyResize) : KeyNone;
    if (c != KeyEscape) return c;

    // A lone Escape and the start of a sequence differ only by what follows,
    // so give the rest of the sequence a short window to arrive.
    int b = readByte(60);
    if (b == KeyNone) return KeyEscape;
    if (b == 'O') return decodeSS3(readByte(60));
    if (b != '[') return KeyEscape;  // Alt+key: not used by the editor

    std::string params;
    int final = KeyNone;
    for (int i = 0; i < 16; ++i) {
      int n = readByte(60);
      if (n == KeyNone) return KeyEscape;
      if (n >= 0x40 && n <= 0x7e) {
        final = n;
        break;
      }
      params.push_back(static_cast<char>(n));
    }
    if (final == KeyNone) return KeyEscape;
    return decodeCSI(params, final);
  }

 private:
  // Reads one byte; timeoutMs < 0 blocks. Returns KeyNone on timeout, and also
  // on EINTR so that a window resize is not mistaken for a key.
  int readByte(int timeoutMs) {
    unsigned char c = 0;
    for (;;) {
      if (timeoutMs >= 0) {
        struct timeval tv;
        tv.tv_sec = timeoutMs / 1000;
        tv.tv_usec = (timeoutMs % 1000) * 1000;
        fd_set set;
        FD_ZERO(&set);
        FD_SET(STDIN_FILENO, &set);
        int ready = select(STDIN_FILENO + 1, &set, nullptr, nullptr, &tv);
        if (ready == 0) return KeyNone;
        if (ready < 0) {
          if (errno == EINTR) return KeyNone;
          return KeyNone;
        }
      }
      ssize_t n = ::read(STDIN_FILENO, &c, 1);
      if (n == 1) return c;
      if (n < 0 && errno == EINTR) return KeyNone;  // most likely SIGWINCH
      if (n == 0) return KeyNone;
    }
  }

  static int decodeSS3(int c) {
    switch (c) {
      case 'A': return KeyArrowUp;
      case 'B': return KeyArrowDown;
      case 'C': return KeyArrowRight;
      case 'D': return KeyArrowLeft;
      case 'H': return KeyHome;
      case 'F': return KeyEnd;
      default: return KeyEscape;
    }
  }

  // params is what sits between "\x1b[" and the final byte: an optional first
  // number, then an optional modifier after a semicolon. The modifier encodes
  // the held keys as 1 + (shift 1, alt 2, ctrl 4), so Shift+Up is "\x1b[1;2A".
  static void splitParams(const std::string& params, int& first, int& modifier) {
    first = 1;
    modifier = 1;
    const size_t semicolon = params.find(';');
    if (semicolon == std::string::npos) {
      if (!params.empty()) first = std::atoi(params.c_str());
    } else {
      if (semicolon > 0) first = std::atoi(params.substr(0, semicolon).c_str());
      const std::string tail = params.substr(semicolon + 1);
      if (!tail.empty()) modifier = std::atoi(tail.c_str());
    }
    if (first <= 0) first = 1;
    if (modifier <= 0) modifier = 1;
  }

  static int decodeCSI(const std::string& params, int final) {
    // rxvt reports the shifted arrows with a lowercase final byte instead of a
    // modifier parameter.
    switch (final) {
      case 'a': return KeyShiftArrowUp;
      case 'b': return KeyShiftArrowDown;
      case 'c': return KeyShiftArrowRight;
      case 'd': return KeyShiftArrowLeft;
      default: break;
    }

    int first = 1;
    int modifier = 1;
    splitParams(params, first, modifier);
    const int held = modifier - 1;
    const bool shift = (held & 1) != 0;
    const bool control = (held & 4) != 0;

    switch (final) {
      case 'A': return control ? KeyCtrlArrowUp : (shift ? KeyShiftArrowUp : KeyArrowUp);
      case 'B': return control ? KeyCtrlArrowDown : (shift ? KeyShiftArrowDown : KeyArrowDown);
      case 'C': return control ? KeyCtrlArrowRight : (shift ? KeyShiftArrowRight : KeyArrowRight);
      case 'D': return control ? KeyCtrlArrowLeft : (shift ? KeyShiftArrowLeft : KeyArrowLeft);
      case 'H': return shift ? KeyShiftHome : KeyHome;
      case 'F': return shift ? KeyShiftEnd : KeyEnd;
      case '~': break;
      default: return KeyEscape;
    }

    switch (first) {
      case 1:
      case 7: return shift ? KeyShiftHome : KeyHome;
      case 3: return KeyDelete;
      case 4:
      case 8: return shift ? KeyShiftEnd : KeyEnd;
      case 5: return shift ? KeyShiftPageUp : KeyPageUp;
      case 6: return shift ? KeyShiftPageDown : KeyPageDown;
      default: return KeyEscape;
    }
  }

  struct termios m_saved {};
  bool m_active = false;
};

// One line of the frame. Text is truncated at the terminal width, counting only
// the characters that actually occupy a column (SGR sequences do not).
class Row {
 public:
  explicit Row(int width) : m_width(width) {}

  // style may combine several SGR sequences, e.g. inverse video over a colour.
  Row& put(const std::string& s, const std::string& style = std::string()) {
    if (m_used >= m_width || s.empty()) return *this;
    const int room = m_width - m_used;
    const std::string clipped = (static_cast<int>(s.size()) <= room) ? s : s.substr(0, static_cast<size_t>(room));
    if (!style.empty()) m_out += style;
    m_out += clipped;
    if (!style.empty()) m_out += sgr::reset;
    m_used += static_cast<int>(clipped.size());
    return *this;
  }

  Row& put(char c, const std::string& style = std::string()) { return put(std::string(1, c), style); }

  // Pads with spaces so that what follows is right-aligned against the edge.
  Row& padTo(int column) {
    if (column > m_used) put(std::string(static_cast<size_t>(column - m_used), ' '));
    return *this;
  }

  int used() const { return m_used; }
  int room() const { return m_width - m_used; }
  const std::string& str() const { return m_out; }

 private:
  std::string m_out;
  int m_width;
  int m_used = 0;
};

// A whole frame: rows are appended in order, then flushed in a single write.
class Frame {
 public:
  explicit Frame(Size size) : m_size(size) { m_out.reserve(8192); }

  void add(const Row& row) {
    if (m_rows >= m_size.rows) return;
    if (m_rows > 0) m_out += "\r\n";
    m_out += row.str();
    m_out += "\x1b[0m\x1b[K";  // reset the style and clear the rest of the line
    m_rows++;
  }

  void addBlank() { add(Row(m_size.cols)); }

  // Fills the frame up, places the caret and pushes everything to the terminal.
  void flush(int cursorRow, int cursorCol, bool showCursor) {
    while (m_rows < m_size.rows) addBlank();
    std::string out = "\x1b[?25l\x1b[H";  // hide the caret while redrawing
    out += m_out;
    char move[32];
    snprintf(move, sizeof(move), "\x1b[%d;%dH", cursorRow + 1, cursorCol + 1);
    out += move;
    if (showCursor) out += "\x1b[?25h";
    Terminal::put(out);
  }

 private:
  Size m_size;
  std::string m_out;
  int m_rows = 0;
};

}  // namespace term

#endif /* end of include guard: TERMINAL_HPP */
