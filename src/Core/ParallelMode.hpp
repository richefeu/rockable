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

#ifndef PARALLEL_MODE_HPP
#define PARALLEL_MODE_HPP

#include <array>
#include <string_view>
#include <optional>

enum ParallelMode {
  DefaultParallelMode,
  InteractionBuffer,  // copy in a intermediate array
  CellMutex,  // one thread = one cell + mutexe
  WaveMethod,  // Use the wave method (27 waves), one thread = one cell, no mutexe
  WaveMethodBlock  // Use the wave method (8 waves), one thread = at least 8 cells, no mutexe
};

inline ParallelMode parallel_mode_from_string(std::string_view s)
{
  if (s == "Default") return DefaultParallelMode;
  if (s == "DefaultParallelMode") return DefaultParallelMode;
  if (s == "InteractionBuffer") return InteractionBuffer;
  if (s == "CellMutex") return CellMutex;
  if (s == "WaveMethod") return WaveMethod;
  if (s == "WaveMethodBlock") return WaveMethodBlock;
  return DefaultParallelMode;
}

static constexpr std::array<std::string_view, 5> names = {
  "DefaultParallelMode",
  "InteractionBuffer",
  "CellMutex",
  "WaveMethod",
  "WaveMethodBlock"
};

inline std::string_view to_string(ParallelMode m)
{
  return names[static_cast<int>(m)];
}

inline std::optional<ParallelMode> from_string(std::string_view s)
{
  for (int i = 0; i < (int)names.size(); ++i) {
      if (names[i] == s)
          return static_cast<ParallelMode>(i);
  }
  return std::nullopt;
}

inline std::istream& operator>>(std::istream& is, ParallelMode& mode)
{
  std::string tmp;
  is >> tmp;  // lit le token texte
  mode = parallel_mode_from_string(tmp);
  return is;
}

#endif
