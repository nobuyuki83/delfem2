/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_DEF_H
#define DFM2_DEF_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/mats.h"

// ---------------------------

namespace delfem2 {

class CDef_SingleLaplacian
{
public:
  void Init(const std::vector<double>& aXYZ0,
            const std::vector<unsigned int>& aTri);
  void Solve(std::vector<double>& aXYZ1,
             const std::vector<double>& aXYZ0,
             const std::vector<int>& aBCFlag);
public:
  CMatrixSparse<double> mat_A;
  std::vector<double> aRhs0, aRhs1;
  std::vector<double> aHistConv;
};

class CDef_LaplacianLinear{
public:
  void Init(const std::vector<double>& aXYZ0,
            const std::vector<unsigned int>& aTri,
            bool is_preconditioner);
  void Solve(std::vector<double>& aXYZ1,
             const std::vector<double>& aXYZ0,
             const std::vector<int>& aBCFlag);
public:
  CMatrixSparse<double> mat_A;
  bool is_preconditioner;
};

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/def.cpp"
#endif

#endif /* def_h */
