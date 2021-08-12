/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_LP_H
#define DFM2_LP_H

#include <math.h>
#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

class CLinPro {
public:
  enum EQ_TYPE {
    LE, GE, EQ
  };
public:
  CLinPro() {}

  void AddEqn(const std::vector<double> &aW, double rhs, EQ_TYPE type);

  int Precomp(int &nitr);

  int Solve(std::vector<double> &solution, double &opt_val, int &nitr,
            const std::vector<double> &aCoeffTrg) const;

  std::vector<double> GetValid() const;

  void Print() const;

private:
  class CEq {
  public:
    bool IsValid(const std::vector<double> &sol, double tol = 1.0e-20) const {
      double sum = -rhs;
      for (unsigned int ic = 0; ic < aCoeff.size(); ++ic) { sum += aCoeff[ic] * sol[ic]; }
      if (itype == EQ) {
        if (fabs(sum) < tol) { return true; }
        else {
          std::cout << "  diff:" << fabs(sum) << std::endl;
          return false;
        }
      } else if (itype == LE) {
        if (sum < tol) { return true; }
        else { return false; }
      } else {
        if (sum > -tol) { return true; }
        else { return false; }
      }
    }

  public:
    std::vector<double> aCoeff;
    double rhs;
    EQ_TYPE itype; // 0:le, 1:ge 2:eq
  };

public:
  std::vector<CEq> aEq;
  //
  unsigned int neq, nvar, nslk, nart;
  std::vector<unsigned int> map_col2row;
  std::vector<double> A;
};

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/lp.cpp"
#endif

#endif
