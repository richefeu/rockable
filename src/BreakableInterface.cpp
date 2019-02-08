// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#include "BreakableInterface.hpp"

BreakableInterface::BreakableInterface()
    : i(0), j(0), kn(0.0), kt(0.0), kr(0.0), fn0(0.0), ft0(0.0), mom0(0.0), power(2.0), dn0(0.0), isInner(1) {}

BreakableInterface::BreakableInterface(size_t I, size_t J)
    : i(I), j(J), kn(0.0), kt(0.0), kr(0.0), fn0(0.0), ft0(0.0), mom0(0.0), power(2.0), dn0(0.0), isInner(1) {}
