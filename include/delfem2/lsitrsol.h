/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @discussion splitting this file into "vecx.h" and "itersol.h" in the future
 */

#ifndef DFM2_LSITRSOL_H
#define DFM2_LSITRSOL_H

#include "delfem2/dfm2_inline.h"
#include <vector>
#include <cassert>
#include <complex>
#include <iostream>

namespace delfem2 {

// ----------------------------------------------------------------------

/**
 * @brief solve linear system using conjugate gradient method
 * @param mat (in)  a template class with member function "MatVec" with  {y} = alpha*[A]{x} + beta*{y}
 */
template<typename REAL, class MAT, class VEC>
std::vector<REAL> Solve_CG(
    VEC& r_vec,
    VEC& u_vec,
    REAL conv_ratio_tol,
    unsigned int max_iteration,
    const MAT& mat,
    VEC& Ap_vec,
    VEC& p_vec)
{
  std::vector<REAL> aConv;
  u_vec.setZero();
  REAL sqnorm_res = r_vec.dot(r_vec);
  if (sqnorm_res < 1.0e-30) { return aConv; }
  REAL inv_sqnorm_res_ini = 1 / sqnorm_res;
  p_vec = r_vec;  // {p} = {r}  (Set Initial Serch Direction)
  for (unsigned int iitr = 0; iitr < max_iteration; iitr++) {
    REAL alpha;
    {  // alpha = (r,r) / (p,Ap)
      AddMatVec(Ap_vec, 0.0, 1.0, mat, p_vec); // {Ap_vec} = [mat]*{p_vec}
      const REAL pAp = p_vec.dot(Ap_vec);
      alpha = sqnorm_res / pAp;
    }
    AddScaledVec(u_vec, alpha, p_vec);    // {x} = +alpha*{ p} + {x} (update x)
    AddScaledVec(r_vec, -alpha, Ap_vec);  // {r} = -alpha*{Ap} + {r}
    const REAL sqnorm_res_new = r_vec.dot(r_vec);
    REAL conv_ratio = sqrt(sqnorm_res_new * inv_sqnorm_res_ini);
    aConv.push_back(conv_ratio);
    if (conv_ratio < conv_ratio_tol) { return aConv; }
    //
    {
      const REAL beta = sqnorm_res_new / sqnorm_res; // beta = (r1,r1) / (r0,r0)
      sqnorm_res = sqnorm_res_new;
      ScaleAndAddVec(p_vec,beta,r_vec); // {p} = {r} + beta*{p}
    }
  }
  return aConv;
}

} // delfem2

  
#endif // DFM2_LSITRSOL
