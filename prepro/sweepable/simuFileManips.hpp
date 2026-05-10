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

#pragma once

/*
===============================================================================
simuFileManips.hpp
-------------------------------------------------------------------------------

A tiny, header-only C++ utility for simple text file manipulation with a
fluent (chainable) API. Designed for simulation workflows such as parameter
sweeps where configuration files must be duplicated and slightly modified.

-------------------------------------------------------------------------------
FEATURES
-------------------------------------------------------------------------------
- Read a text file into memory
- Find and replace lines starting with a prefix
- Replace multiple matching lines
- Replace tokens inside a line
- Insert and append lines
- Create folders and save files
- Custom output filename
- Optional "silent" mode for missing matches
- Fully chainable API (fluent interface)

-------------------------------------------------------------------------------
BASIC USAGE
-------------------------------------------------------------------------------

SFManip()
    .createFolder("dt%g", dt)
    .read("input.txt")
    .findLineStartingWith("dt ")
    .replaceBy("dt %f", dt)
    .saveInFolder();

-------------------------------------------------------------------------------
METHOD OVERVIEW
-------------------------------------------------------------------------------

read(filename)
    Load a file into memory.

findLineStartingWith(prefix)
    Select the first matching line.

replaceBy(fmt, ...)
    Replace selected line.

replaceAllStartingWith(prefix, fmt, ...)
    Replace all matching lines.

replaceInLine(target, replacement)
    Replace substring inside selected line.

insertAfter(prefix, newLine)
    Insert a line after a match.

appendLine(line)
    Append a line at the end.

createFolder(fmt, ...)
    Create and store output folder.

saveInFolder()
saveInFolder(fmt, ...)
    Save file to disk.

setOutputFilename(name)
    Change output filename (default: input.txt).

silentIfNotFound(bool)
    Disable exceptions for missing matches.

-------------------------------------------------------------------------------
NOTES
-------------------------------------------------------------------------------
- Uses printf-style formatting
- Throws std::runtime_error unless silent mode is enabled
- Requires C++17 (<filesystem>)

===============================================================================
*/

#include <cstdarg>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

class SFManip {
 private:
  std::vector<std::string> lines;
  int currentLineIndex = -1;

  std::string currentFolder;
  std::string outputFilename = "input.txt";
  bool ignoreMissing = false;

  // --- formatting helper ---
  std::string vformat(const char* fmt, va_list args) {
    va_list tmp;
    va_copy(tmp, args);
    int size = std::vsnprintf(nullptr, 0, fmt, tmp);
    va_end(tmp);

    std::vector<char> buf(size + 1);
    std::vsnprintf(buf.data(), buf.size(), fmt, args);
    return std::string(buf.data());
  }

 public:
  SFManip() = default;
  bool verbose{true};

  // ------------------------------------------------------------------------
  // Collection
  // ------------------------------------------------------------------------
  inline static std::vector<std::string> Collection;

  static void clearCollection() { Collection.clear(); }

  static void saveCollection(const char* path) {

    std::ofstream out(path);

    if (!out) {
      throw std::runtime_error(std::string("Cannot write collection file: ") + path);
    }

    for (const auto& p : Collection) {
      out << p << "\n";
    }
  }

  // ------------------------------------------------------------------------
  // File I/O
  // ------------------------------------------------------------------------

  // Load file into memory
  SFManip& read(const std::string& filename) {
    lines.clear();
    currentLineIndex = -1;

    std::ifstream file(filename);
    if (!file) {
      throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
      lines.push_back(line);
    }

    if (verbose) {
      std::cout << " > read " << filename << std::endl;
    }

    return *this;
  }

  // Set output filename (default: input.txt)
  SFManip& setOutputFilename(const std::string& name) {
    outputFilename = name;

    if (verbose) {
      std::cout << " > setOutputFilename " << name << std::endl;
    }

    return *this;
  }

  // ------------------------------------------------------------------------
  // Line selection / modification
  // ------------------------------------------------------------------------

  // Find first line starting with prefix
  SFManip& findLineStartingWith(const std::string& prefix) {
    currentLineIndex = -1;

    for (size_t i = 0; i < lines.size(); ++i) {
      if (lines[i].rfind(prefix, 0) == 0) {
        currentLineIndex = static_cast<int>(i);
        break;
      }
    }

    if (currentLineIndex == -1 && !ignoreMissing) {
      throw std::runtime_error("Line starting with '" + prefix + "' not found.");
    }

    if (verbose) {
      std::cout << " > findLineStartingWith " << prefix << std::endl;
    }

    return *this;
  }

  // Replace selected line
  SFManip& replaceBy(const char* fmt, ...) {
    if (currentLineIndex == -1) {
      if (!ignoreMissing) throw std::runtime_error("No line selected.");
      return *this;
    }

    va_list args;
    va_start(args, fmt);
    lines[currentLineIndex] = vformat(fmt, args);
    va_end(args);

    if (verbose) {
      std::cout << " > replaceBy " << lines[currentLineIndex] << std::endl;
    }

    return *this;
  }

  // Replace all lines starting with prefix
  SFManip& replaceAllStartingWith(const std::string& prefix, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::string replacement = vformat(fmt, args);
    va_end(args);

    bool found = false;

    for (auto& line : lines) {
      if (line.rfind(prefix, 0) == 0) {
        line = replacement;
        found = true;
      }
    }

    if (!found && !ignoreMissing) {
      throw std::runtime_error("No lines starting with '" + prefix + "' found.");
    }

    if (verbose) {
      std::cout << " > replaceAllStartingWith " << prefix << "  --> " << replacement << std::endl;
    }

    return *this;
  }

  // Replace substring inside selected line
  SFManip& replaceInLine(const std::string& target, const std::string& replacement) {
    if (currentLineIndex == -1) {
      if (!ignoreMissing) throw std::runtime_error("No line selected.");
      return *this;
    }

    auto& line = lines[currentLineIndex];
    size_t pos = line.find(target);

    if (pos == std::string::npos) {
      if (!ignoreMissing) throw std::runtime_error("Target not found in line.");
      return *this;
    }

    line.replace(pos, target.length(), replacement);
    
    if (verbose) {
      std::cout << " > replaceInLine: " << target << " --> " << replacement << std::endl;
    }
    
    return *this;
  }

  // Insert a line after a match
  SFManip& insertAfter(const std::string& prefix, const std::string& newLine) {
    if (verbose) {
      std::cout << " > insertAfter: " << prefix << "  -->  " << newLine << std::endl;
    }
    
    for (size_t i = 0; i < lines.size(); ++i) {
      if (lines[i].rfind(prefix, 0) == 0) {
        lines.insert(lines.begin() + i + 1, newLine);
        return *this;
      }
    }

    if (!ignoreMissing) {
      throw std::runtime_error("Line starting with '" + prefix + "' not found.");
    }

    return *this;
  }

  // Append a line at end of file
  SFManip& appendLine(const std::string& line) {
    lines.push_back(line);
    if (verbose) {
      std::cout << " > appendLine: " << line << std::endl;
    }
    
    return *this;
  }

  // ------------------------------------------------------------------------
  // Folder / output handling
  // ------------------------------------------------------------------------

  // Create folder and store it
  SFManip& createFolder(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    currentFolder = vformat(fmt, args);
    va_end(args);

    std::filesystem::create_directories(currentFolder);
    if (verbose) {
      std::cout << " > create_directories: " << currentFolder << std::endl;
    }
    return *this;
  }

  // Save using stored folder
  SFManip& saveInFolder() {
    if (currentFolder.empty()) {
      throw std::runtime_error("No folder defined.");
    }

    std::string outputPath = currentFolder + "/" + outputFilename;
    std::ofstream out(outputPath);

    if (!out) {
      throw std::runtime_error("Cannot write file: " + outputPath);
    }

    for (const auto& line : lines) {
      out << line << "\n";
    }

    if (verbose) {
      std::cout << " > saveInFolder " << outputPath << std::endl;
    }

    Collection.push_back(outputPath);

    return *this;
  }

  // Stateless save
  SFManip& saveInFolder(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::string folder = vformat(fmt, args);
    va_end(args);

    std::filesystem::create_directories(folder);

    std::string outputPath = folder + "/" + outputFilename;
    std::ofstream out(outputPath);

    if (!out) {
      throw std::runtime_error("Cannot write file: " + outputPath);
    }

    for (const auto& line : lines) {
      out << line << "\n";
    }

    if (verbose) {
      std::cout << " > saveInFolder " << outputPath << std::endl;
    }

    Collection.push_back(outputPath);

    return *this;
  }

  // ------------------------------------------------------------------------
  // Behavior control
  // ------------------------------------------------------------------------

  // Enable/disable silent mode (ignore missing matches)
  SFManip& silentIfNotFound(bool value = true) {
    ignoreMissing = value;
    if (verbose) {
      std::cout << " > silentIfNotFound" << std::endl;
    }
    return *this;
  }
};