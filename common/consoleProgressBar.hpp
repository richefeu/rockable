#ifndef CONSOLE_PROGRESS_BAR
#define CONSOLE_PROGRESS_BAR

#include <iomanip>
#include <iostream>
#include <string>

class ConsoleProgressBar {
private:
  size_t nmax;
  size_t width;
  std::string title;
  char progressChar;
  char voidChar;
  char openChar;
  char closeChar;

  char validChar(char c, char cdef) {
    int ci = (int)c;
    if (ci >= 33 && ci <= 127) {
      return (char)c;
    }
    return (char)cdef;
  }

public:
  ConsoleProgressBar(size_t n = 100) {
    nmax = n;
    width = 50;
    title = "";
    progressChar = '|';
    voidChar = ' ';
    openChar = '[';
    closeChar = ']';
  }

  void setMax(size_t n) { nmax = n; }
  void setTitle(const char *t) { title = t; }
  void setWidth(size_t w) { width = w; }

  void setProgressChar(char c) { progressChar = (char)validChar(c, '|'); }
  void setVoidChar(char c) { voidChar = (char)validChar(c, ' '); }
  void setOpenChar(char c) { openChar = (char)validChar(c, '['); }
  void setCloseCharChar(char c) { closeChar = (char)validChar(c, ']'); }

  // Process has done i out of n rounds,
  // and we want a bar of width w and resolution r.
  void update(size_t x, std::ostream &os = std::cerr) {
    if ((x != nmax) && (x % (nmax / 100 + 1) != 0))
      return;

    float ratio = x / (float)nmax;
    size_t c = ratio * width;

    os << title << std::setw(3) << (size_t)(ratio * 100) << "% " << openChar;
    for (size_t x = 0; x < c; x++) {
      os << progressChar;
    }
    for (size_t x = c; x < width; x++) {
      os << voidChar;
    }
    os << closeChar << '\r' << std::flush;
  }
};

#if 0
int main(int argc, char const *argv[]) {
  ConsoleProgressBar cpb(100);
  cpb.setTitle("Title: ");
  cpb.setWidth(50);
  cpb.setProgressChar(124);
  cpb.setVoidChar('-');
  cpb.update(10, std::cout);
  return 0;
}
#endif

#endif /* end of include guard: CONSOLE_PROGRESS_BAR */
