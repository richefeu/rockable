//        Rockable, 3D-DEM with sphero-polyhedra
//        Copyright (C) 2016-2019  <vincent.richefeu@3sr-grenoble.fr>
//        
//        This program is free software: you can redistribute it and/or modify
//        it under the terms of the GNU General Public License as published by
//        the Free Software Foundation, either version 3 of the License, or
//        (at your option) any later version.
//        
//        This program is distributed in the hope that it will be useful,
//        but WITHOUT ANY WARRANTY; without even the implied warranty of
//        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//        GNU General Public License for more details.
//        
//        You should have received a copy of the GNU General Public License
//        along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "BreakableInterface.hpp"

BreakableInterface::BreakableInterface()
    : i(0), j(0), kn(0.0), kt(0.0), kr(0.0), fn0(0.0), ft0(0.0), mom0(0.0), power(2.0), dn0(0.0), isInner(1) {}

BreakableInterface::BreakableInterface(size_t I, size_t J)
    : i(I), j(J), kn(0.0), kt(0.0), kr(0.0), fn0(0.0), ft0(0.0), mom0(0.0), power(2.0), dn0(0.0), isInner(1) {}
