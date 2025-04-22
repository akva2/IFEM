//==============================================================================
//!
//! \file ASMuSquareGo.h
//!
//! \date Oct 2 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Default bi-unit square LR patch for unit-testing.
//!
//==============================================================================

#ifndef ASM_USQUARE_GO_H
#define ASM_USQUARE_GO_H

#include "GoLR/ASMu2DGo.h"
#include <sstream>

class ASMuSquareGo : public ASMu2DGo
{
public:
  static constexpr const char* square =
    "# LRSPLINE SURFACE\n"
    "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n"
    "\t2\t2\t4\t4\t1\t2\t0\n"
    "# Basis functions:\n"
    "0: [0 0 1 ] x [0 0 1 ] 0 0 (1)\n"
    "1: [0 1 1 ] x [0 0 1 ] 1 0 (1)\n"
    "2: [0 0 1 ] x [0 1 1 ] 0 1 (1)\n"
    "3: [0 1 1 ] x [0 1 1 ] 1 1 (1)\n"
    "# Mesh lines:\n0 x [0, 1] (2)\n"
    "1 x [0, 1] (2)\n"
    "[0, 1] x 0 (2)\n"
    "[0, 1] x 1 (2)\n"
    "# Elements:\n"
    "0 [2] : (0, 0) x (1, 1)    {0, 1, 2, 3}\n";

  explicit ASMuSquareGo(unsigned char n_f = 2,
                        double xshift = 0.0,
                        double yshift = 0.0)
    : ASMu2DGo(2,n_f)
  {
    auto&& addLine = [](double x, double y)
    {
      return std::to_string(x) + " " + std::to_string(y) + '\n';
    };
    std::stringstream geo("200 1 0 0\n2 0\n2 2\n0 0 1 1\n2 2\n0 0 1 1\n" +
                          addLine(xshift,     yshift) +
                          addLine(xshift+1.0, yshift) +
                          addLine(xshift,     yshift+1.0) +
                          addLine(xshift+1.0, yshift+1.0));
    this->read(geo);
  }
  virtual ~ASMuSquareGo() {}
};

#endif
