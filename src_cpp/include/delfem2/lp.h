/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef LP_H
#define LP_H

#include <vector>



class CLinPro
{
public:
  enum EQ_TYPE{ LE, GE, EQ };
public:
  CLinPro(){}
  void AddEqn(const std::vector<double>& aW, double rhs, EQ_TYPE type);
  int Precomp(int& nitr);
  bool Solve(std::vector<double>& solution, double& opt_val, int& nitr,
             const std::vector<double>& aCoeffTrg) const;
//  void SetTarget(const std::vector<double>& aW);
  ////
//  bool IsViable(const std::vector<double>& aW) const;
//  std::vector<double> GetViableSolution(std::vector<double>& sol, int& nitr) const;
  ////
private:
  class CEq
  {
  public:
    std::vector<double> aCoeff;
    double rhs;
    EQ_TYPE itype; // 0:le, 1:ge 2:eq
  };
public:
  std::vector<CEq> aEq;
  ////
  unsigned int neq, nvar, nslk, nart;
  std::vector<int> map_col2row;
  std::vector<double> A;
};


#endif
