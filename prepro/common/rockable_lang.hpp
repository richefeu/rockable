// rockable_lang.hpp -- loads rockable.lang, the shared description of the
// Rockable input-file language (keywords, types, documentation, snippets).
//
// Header-only and dependent on the standard library only, so that any of the
// prepro tools can use it. The tables it fills are the ones that used to be
// hard-coded in confedit.cpp; keeping them in a data file means a new keyword
// can be documented without recompiling anything.
//
// Typical use:
//   rockable_lang::Language lang;
//   lang.load(rockable_lang::findFile(argv[0]));
//   if (lang.kindOf("tmax") == rockable_lang::Kind::Keyword) ...

#ifndef ROCKABLE_LANG_HPP
#define ROCKABLE_LANG_HPP

#include <cstdlib>
#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace rockable_lang {

enum class Kind {
  None,     // unknown word, no highlighting
  Keyword,  // a command word
  Type,     // an enumerated value
  Doc       // documented, but not highlighted
};

struct Snippet {
  std::string title;
  std::string body;
};

class Language {
 public:
  // Parses the given file. Returns false (leaving the tables empty) if it
  // cannot be opened; the caller decides whether that is fatal.
  bool load(const std::string& filename) {
    std::ifstream file(filename.c_str());
    if (!file) return false;

    m_path = filename;
    m_kinds.clear();
    m_docs.clear();
    m_snippets.clear();

    std::string line;
    std::string currentWord;     // entry the doc lines attach to
    std::string pendingDoc;      // doc collected so far
    std::string pendingBlanks;   // blank lines held back until more doc arrives
    bool inDoc = false;
    bool inSnippet = false;
    std::string snippetTitle;
    std::string snippetBody;

    auto flushDoc = [&]() {
      if (inDoc && !currentWord.empty() && !pendingDoc.empty()) {
        m_docs[currentWord] = pendingDoc;
      }
      inDoc = false;
      pendingDoc.clear();
      pendingBlanks.clear();
    };

    auto flushSnippet = [&]() {
      if (inSnippet) {
        while (!snippetBody.empty() && snippetBody.back() == '\n') snippetBody.pop_back();
        if (!snippetBody.empty()) m_snippets.push_back({snippetTitle, snippetBody});
      }
      inSnippet = false;
      snippetTitle.clear();
      snippetBody.clear();
    };

    while (std::getline(file, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();

      // A tag at column 0 always closes whatever block was open.
      if (!line.empty() && line[0] == '[') {
        flushDoc();
        flushSnippet();
        size_t close = line.find(']');
        if (close == std::string::npos) continue;
        std::string tag = line.substr(1, close - 1);
        std::string rest = trim(line.substr(close + 1));
        if (rest.empty()) continue;

        if (tag == "snippet") {
          inSnippet = true;
          snippetTitle = rest;
          currentWord.clear();
        } else {
          currentWord = rest;
          if (tag == "keyword") {
            m_kinds[rest] = Kind::Keyword;
          } else if (tag == "type") {
            m_kinds[rest] = Kind::Type;
          } else if (tag == "doc") {
            m_kinds[rest] = Kind::Doc;
          }
        }
        continue;
      }

      if (inSnippet) {
        snippetBody += line;
        snippetBody += '\n';
        continue;
      }

      // "doc:" opens a documentation block; its first line may follow on the
      // same line. Continuation lines are indented by exactly two spaces.
      if (!inDoc && line.compare(0, 4, "doc:") == 0) {
        inDoc = true;
        pendingDoc = trim(line.substr(4));
        pendingBlanks.clear();
        continue;
      }

      if (inDoc) {
        if (line.empty()) {
          pendingBlanks += '\n';  // kept only if more doc follows
          continue;
        }
        if (line.compare(0, 2, "  ") == 0) {
          if (!pendingDoc.empty()) pendingDoc += '\n';
          pendingDoc += pendingBlanks;
          pendingBlanks.clear();
          pendingDoc += line.substr(2);
          continue;
        }
        flushDoc();
        // fall through: this line is ordinary content again
      }

      // Outside any block, '#' introduces a comment and blanks are ignored.
    }
    flushDoc();
    flushSnippet();

    m_loaded = true;
    return true;
  }

  bool loaded() const { return m_loaded; }
  const std::string& path() const { return m_path; }

  Kind kindOf(const std::string& word) const {
    auto it = m_kinds.find(word);
    return (it == m_kinds.end()) ? Kind::None : it->second;
  }

  // Returns nullptr when the word is not documented.
  const std::string* doc(const std::string& word) const {
    auto it = m_docs.find(word);
    return (it == m_docs.end()) ? nullptr : &it->second;
  }

  const std::vector<Snippet>& snippets() const { return m_snippets; }
  size_t wordCount() const { return m_kinds.size(); }
  size_t docCount() const { return m_docs.size(); }

 private:
  static std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return std::string();
    size_t e = s.find_last_not_of(" \t");
    return s.substr(b, e - b + 1);
  }

  std::unordered_map<std::string, Kind> m_kinds;
  std::unordered_map<std::string, std::string> m_docs;
  std::vector<Snippet> m_snippets;
  std::string m_path;
  bool m_loaded = false;
};

// Looks for rockable.lang in, by order of decreasing priority: the ROCKABLE_LANG
// environment variable, the current directory, the directory holding the running
// executable and its three parents, and finally ~/.rockable. Each of those
// directories is tried both directly and through a "common" subdirectory, which
// is where the file sits in the source tree (prepro/common). Returns an empty
// string when nothing is found.
inline std::string findFile(const char* argv0 = nullptr) {
  auto exists = [](const std::string& p) {
    std::ifstream f(p.c_str());
    return static_cast<bool>(f);
  };

  if (const char* env = std::getenv("ROCKABLE_LANG")) {
    if (exists(env)) return env;
  }

  std::vector<std::string> directories;
  directories.push_back(".");

  if (argv0 != nullptr) {
    const std::string exe(argv0);
    const size_t slash = exe.find_last_of('/');
    const std::string dir = (slash == std::string::npos) ? std::string(".") : exe.substr(0, slash);
    directories.push_back(dir);
    directories.push_back(dir + "/..");
    directories.push_back(dir + "/../..");
    directories.push_back(dir + "/../../..");
  }

  if (const char* home = std::getenv("HOME")) {
    directories.push_back(std::string(home) + "/.rockable");
  }

  for (const std::string& directory : directories) {
    std::string candidate = directory + "/rockable.lang";
    if (exists(candidate)) return candidate;
    candidate = directory + "/common/rockable.lang";  // the prepro/common layout
    if (exists(candidate)) return candidate;
  }
  return std::string();
}

}  // namespace rockable_lang

#endif /* end of include guard: ROCKABLE_LANG_HPP */
