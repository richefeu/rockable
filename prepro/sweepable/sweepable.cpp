//  Copyright or © or Copr. Rockable
//
//  vincent.richefeu@3sr-grenoble.fr
//
//  This software is a computer program whose purpose is
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces.
//  It is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

// sweepable.cpp - Command processor for SFManip

#include "simuFileManips.hpp"

#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

// ============================================================
// Value
// Holds either a double or a string
// ============================================================

struct Value {
  std::variant<double, std::string> v;

  Value() : v("") {}
  Value(double d) : v(d) {}
  Value(const std::string& s) : v(s) {}

  bool isDouble() const { return std::holds_alternative<double>(v); }

  bool isString() const { return std::holds_alternative<std::string>(v); }

  double asDouble() const { return std::get<double>(v); }

  std::string asString() const { return std::get<std::string>(v); }

  std::string toString() const {
    if (isDouble()) {
      char buf[64];
      snprintf(buf, sizeof(buf), "%g", asDouble());
      return buf;
    }
    return asString();
  }
};

// ============================================================
// Context
// Stores variables only
// Loop state is now handled by recursion in CmdProcessor
// ============================================================

struct Context {
  std::map<std::string, Value> vars;

  void set(const std::string& name, const Value& val) { vars[name] = val; }

  bool has(const std::string& name) const { return vars.count(name) > 0; }

  Value get(const std::string& name) const {
    auto it = vars.find(name);

    if (it == vars.end()) {
      throw std::runtime_error("Undefined variable: " + name);
    }

    return it->second;
  }
};

// ============================================================
// Tokenizer
// Splits a line into tokens while respecting quotes
// ============================================================

std::vector<std::string> tokenize(const std::string& line) {
  std::vector<std::string> tokens;

  std::string token;
  bool inQuotes = false;

  for (size_t i = 0; i < line.size(); i++) {
    char c = line[i];

    // Handle escaped quotes
    if (c == '\\' && i + 1 < line.size()) {
      token += line[i + 1];
      i++;
      continue;
    }

    if (c == '"') {
      inQuotes = !inQuotes;
      continue;
    }

    if (isspace(static_cast<unsigned char>(c)) && !inQuotes) {
      if (!token.empty()) {
        tokens.push_back(token);
        token.clear();
      }
    } else {
      token += c;
    }
  }

  if (!token.empty()) {
    tokens.push_back(token);
  }

  return tokens;
}

// ============================================================
// Variable Resolution
// Replaces:
//   $var or ${var}
// with variable contents
// ============================================================

std::string resolve(const std::string& arg, const Context& ctx) {
  if (arg.empty()) {
    return arg;
  }

  std::string result;
  size_t pos = 0;

  while (pos < arg.size()) {

    size_t dollar = arg.find('$', pos);

    if (dollar == std::string::npos) {
      result += arg.substr(pos);
      break;
    }

    // Copy text before $
    result += arg.substr(pos, dollar - pos);

    size_t start = dollar + 1;

    if (start >= arg.size()) {
      result += "$";
      break;
    }

    std::string name;

    // ----------------------------------------------------
    // ${var} syntax
    // ----------------------------------------------------
    if (arg[start] == '{') {

      size_t endBrace = arg.find('}', start);

      if (endBrace == std::string::npos) {
        // malformed, treat literally
        result += arg.substr(dollar, 1);
        pos = start;
        continue;
      }

      name = arg.substr(start + 1, endBrace - start - 1);
      pos = endBrace + 1;
    }

    // ----------------------------------------------------
    // $var syntax
    // ----------------------------------------------------
    else {

      size_t end = start;

      while (end < arg.size() && (isalnum(static_cast<unsigned char>(arg[end])) || arg[end] == '_')) {
        end++;
      }

      name = arg.substr(start, end - start);
      pos = end;
    }

    // ----------------------------------------------------
    // Replace variable
    // ----------------------------------------------------
    if (ctx.has(name)) {
      result += ctx.get(name).toString();
    } else {
      // keep original text if undefined
      result += arg.substr(dollar, pos - dollar);
    }
  }

  return result;
}

// ============================================================
// Command Processor
// ============================================================

class CmdProcessor {
 private:
  Context ctx;
  SFManip sf;

  std::vector<std::string> commandFileLines;

  // Program counter
  size_t pc = 0;

 private:
  // --------------------------------------------------------
  // Error helper
  // --------------------------------------------------------

  [[noreturn]]
  void error(const std::string& msg) {
    throw std::runtime_error(msg);
  }

  // --------------------------------------------------------
  // Parse string into Value
  // --------------------------------------------------------

  Value parseValue(const std::string& s) {
    try {
      size_t pos = 0;

      double d = stod(s, &pos);

      // Entire string parsed as number
      if (pos == s.size()) {
        return Value(d);
      }

    } catch (...) {
    }

    return Value(s);
  }

  // --------------------------------------------------------
  // Resolve all variables in argument list
  // --------------------------------------------------------

  void resolveArgs(std::vector<std::string>& args) {
    for (auto& a : args) {
      a = resolve(a, ctx);
    }
  }

  // --------------------------------------------------------
  // Execute commands in range [start, end)
  // --------------------------------------------------------

  void executeRange(size_t start, size_t end) {
    pc = start;

    while (pc < end) {

      std::vector<std::string> tokens = tokenize(commandFileLines[pc]);

      pc++;

      if (tokens.empty()) {
        continue;
      }

      const std::string& cmd = tokens[0];

      // -------------------------
      // foreach
      // -------------------------

      if (cmd == ".foreach") {
        processForeach(tokens);
      }

      // -------------------------
      // endforeach
      // -------------------------

      else if (cmd == ".endforeach") {
        // handled by processForeach()
        return;
      }

      // -------------------------
      // set
      // -------------------------

      else if (cmd == ".set") {
        processSet(tokens);
      }

      // -------------------------
      // verbosity
      // -------------------------

      else if (cmd == ".verbose") {
        sf.verbose = true;
      } else if (cmd == ".mute" || cmd == ".quiet") {
        sf.verbose = false;
      }

      // -------------------------
      // SFManip command
      // -------------------------

      else if (!cmd.empty() && cmd[0] == '.') {
        processSFCommand(cmd.substr(1), tokens);
      }
    }
  }

  // --------------------------------------------------------
  // Find matching .endforeach
  // Supports nested loops correctly
  // --------------------------------------------------------

  size_t findMatchingEndforeach(size_t foreachPc) {

    int depth = 1;

    for (size_t i = foreachPc + 1; i < commandFileLines.size(); i++) {

      auto t = tokenize(commandFileLines[i]);

      if (t.empty()) {
        continue;
      }

      if (t[0] == ".foreach") {
        depth++;
      } else if (t[0] == ".endforeach") {
        depth--;

        if (depth == 0) {
          return i;
        }
      }
    }

    error("Unmatched .foreach");

    return 0;
  }

  // --------------------------------------------------------
  // Process .foreach
  // Fully recursive
  // Proper nested loop support
  // --------------------------------------------------------

  void processForeach(const std::vector<std::string>& tokens) {

    if (tokens.size() < 3) {
      error(".foreach requires variable and values");
    }

    std::string var = tokens[1];

    // Parse values
    std::vector<Value> values;

    for (size_t i = 2; i < tokens.size(); i++) {
      values.push_back(parseValue(tokens[i]));
    }

    // Current .foreach line index
    size_t foreachLine = pc - 1;

    // Find matching .endforeach
    size_t endforeachLine = findMatchingEndforeach(foreachLine);

    // Loop body range
    size_t bodyStart = foreachLine + 1;
    size_t bodyEnd = endforeachLine;

    // Save old variable value if it exists
    bool hadOldValue = ctx.has(var);
    Value oldValue;

    if (hadOldValue) {
      oldValue = ctx.get(var);
    }

    // Execute loop
    for (const auto& v : values) {

      ctx.set(var, v);

      executeRange(bodyStart, bodyEnd);
    }

    // Restore previous variable
    if (hadOldValue) {
      ctx.set(var, oldValue);
    }

    // Continue after endforeach
    pc = endforeachLine + 1;
  }

  // --------------------------------------------------------
  // Process .set
  // --------------------------------------------------------

  void processSet(const std::vector<std::string>& tokens) {

    if (tokens.size() < 3) {
      error(".set requires name and value");
    }

    std::string name = tokens[1];

    Value val = parseValue(tokens[2]);

    ctx.set(name, val);
  }

  // --------------------------------------------------------
  // Process SFManip commands
  // --------------------------------------------------------

  void processSFCommand(const std::string& cmd, const std::vector<std::string>& tokens) {
    std::vector<std::string> args(tokens.begin() + 1, tokens.end());

    resolveArgs(args);

    // ----------------------------------------------------
    // read
    // ----------------------------------------------------

    if (cmd == "read") {

      if (args.size() < 1) {
        error("read requires filename");
      }

      sf.read(args[0]);
    }

    // ----------------------------------------------------
    // setOutputFilename
    // ----------------------------------------------------

    else if (cmd == "setOutputFilename") {

      if (args.size() < 1) {
        error("setOutputFilename requires name");
      }

      sf.setOutputFilename(args[0]);
    }

    // ----------------------------------------------------
    // findLineStartingWith
    // ----------------------------------------------------

    else if (cmd == "findLineStartingWith") {

      if (args.size() < 1) {
        error("findLineStartingWith requires prefix");
      }

      sf.findLineStartingWith(args[0]);
    }

    // ----------------------------------------------------
    // replaceBy
    // ----------------------------------------------------

    else if (cmd == "replaceBy") {

      if (args.size() < 1) {
        error("replaceBy requires format");
      }

      sf.replaceBy(args[0].c_str());
    }

    // ----------------------------------------------------
    // replaceAllStartingWith
    // ----------------------------------------------------

    else if (cmd == "replaceAllStartingWith") {
      if (args.size() < 2) {
        error(
            "replaceAllStartingWith requires "
            "prefix and format");
      }
      sf.replaceAllStartingWith(args[0], args[1].c_str());
    }

    // ----------------------------------------------------
    // replaceInLine
    // ----------------------------------------------------

    else if (cmd == "replaceInLine") {
      if (args.size() < 2) {
        error(
            "replaceInLine requires "
            "target and replacement");
      }
      sf.replaceInLine(args[0], args[1]);
    }

    // ----------------------------------------------------
    // insertAfter
    // ----------------------------------------------------

    else if (cmd == "insertAfter") {
      if (args.size() < 2) {
        error(
            "insertAfter requires "
            "prefix and line");
      }
      sf.insertAfter(args[0], args[1]);
    }

    // ----------------------------------------------------
    // appendLine
    // ----------------------------------------------------

    else if (cmd == "appendLine") {
      if (args.size() < 1) {
        error("appendLine requires line");
      }
      sf.appendLine(args[0]);
    }

    // ----------------------------------------------------
    // createFolder
    // ----------------------------------------------------

    else if (cmd == "createFolder") {
      if (args.size() < 1) {
        error("createFolder requires format");
      }
      sf.createFolder(args[0].c_str());
    }

    // ----------------------------------------------------
    // saveInFolder
    // ----------------------------------------------------

    else if (cmd == "saveInFolder") {
      if (args.empty()) {
        sf.saveInFolder();
      } else {
        sf.saveInFolder(args[0].c_str());
      }
    }

    // ----------------------------------------------------
    // silentIfNotFound
    // ----------------------------------------------------

    else if (cmd == "silentIfNotFound") {
      bool val = args.empty() || args[0] != "false";

      sf.silentIfNotFound(val);
    }

    // ----------------------------------------------------
    // Collection
    // ----------------------------------------------------

    else if (cmd == "saveCollection") {
      SFManip::saveCollection(args[0].c_str());
      if (sf.verbose == true) {
        std::cout << " > saveCollection " << args[0] << std::endl;
      }
    } else if (cmd == "clearCollection") {
      SFManip::clearCollection();
      if (sf.verbose == true) {
        std::cout << " > clearCollection" << std::endl;
      }
    }

    // ----------------------------------------------------
    // unknown command
    // ----------------------------------------------------

    else {
      error("Unknown command: " + cmd);
    }
  }

 public:
  // --------------------------------------------------------
  // Load command file
  // --------------------------------------------------------

  void loadCommands(const std::string& filename) {

    std::ifstream f(filename);

    if (!f) {
      error("Cannot open command file: " + filename);
    }

    std::string line;

    while (getline(f, line)) {

      std::string cleaned;

      bool inQuotes = false;

      for (size_t i = 0; i < line.size(); i++) {

        char c = line[i];

        // Handle escaped chars
        if (c == '\\' && i + 1 < line.size()) {
          cleaned += c;
          cleaned += line[i + 1];
          i++;
          continue;
        }

        if (c == '"') {
          inQuotes = !inQuotes;
        }

        // Strip comments outside quotes
        if (c == '#' && !inQuotes) {
          break;
        }

        cleaned += c;
      }

      // Trim trailing whitespace
      while (!cleaned.empty() && isspace(static_cast<unsigned char>(cleaned.back()))) {
        cleaned.pop_back();
      }

      if (!cleaned.empty()) {
        commandFileLines.push_back(cleaned);
      }
    }
  }

  // --------------------------------------------------------
  // Run script
  // --------------------------------------------------------

  void run() { executeRange(0, commandFileLines.size()); }
};

// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[]) {

  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <command_file>" << std::endl;

    std::ofstream file("batchable.sh");

    // clang-format off
    file   
    << "#!/bin/sh\n"
    << "\n"
    << "INPUT_LIST=\"$1\"\n"
    << "COMMAND=\"$2\"\n"
    << "NPROCS=\"$3\"\n"
    << "\n"
    << "if [ $# -ne 3 ]; then\n"
    << "    echo \"Usage: $0 input_list command nprocs\"\n"
    << "    exit 1\n"
    << "fi\n"
    << "\n"
    << "xargs -I{} -P \"$NPROCS\" sh -c '\n"
    << "    input=\"$1\"\n"
    << "    command=\"$2\"\n"
    << "\n"
    << "    dir=$(dirname \"$input\")\n"
    << "    file=$(basename \"$input\")\n"
    << "\n"
    << "    echo \"Running $file in folder $dir\"\n"
    << "\n"
    << "    cd \"$dir\" || {\n"
    << "        echo \"ERROR: cannot cd into $dir\" >&2\n"
    << "        exit 1\n"
    << "    }\n"
    << "\n"
    << "    eval \"$command \\\"$file\\\"\" > run.log 2>&1\n"
    << "' _ {} \"$COMMAND\" < \"$INPUT_LIST\"\n";

    // clang-format on

    std::cout << "File 'batchable.sh' has been created" << std::endl;

    return 1;
  }

  try {

    CmdProcessor proc;
    proc.loadCommands(argv[1]);
    proc.run();

  } catch (const std::exception& e) {

    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}