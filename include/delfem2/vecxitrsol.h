/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @discussion splitting this file into "vecx.h" and "itersol.h" in the future
 */

#ifndef DFM2_VECXITRSOL_H
#define DFM2_VECXITRSOL_H

#include <vector>
#include <cassert>
#include <complex>
#include <iostream>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 *@brief inner prodcut of a vector
 *@details defiend for "double", "float", "std::complex<double>". Inner product of complex value using conjugate value
 */
template<typename T>
DFM2_INLINE T DotX(
    const T *va,
    const T *vb,
    size_t n);

template<typename T>
DFM2_INLINE T Dot(
    const std::vector<T> &r_vec,
    const std::vector<T> &u_vec);

/**
 * @brief Eucledian distance
 */
template<typename T>
T Distance(
    const std::vector<T> &r_vec,
    const std::vector<T> &u_vec);


template <typename T>
void ScaleX(
    T *p0,
    T s,
    size_t n);

template<typename T>
void XPlusAY(
    std::vector<T> &X,
    size_t nDoF,
    const std::vector<int> &aBCFlag,
    T alpha,
    const std::vector<T> &Y);

template<typename T>
void AXPY(
    T a,
    const std::vector<T> &x,
    std::vector<T> &y);

template<typename T>
void AXPY(
    T a,
    const T *x,
    T *y,
    unsigned int n);

/**
 * @brief set zero to the element if the flag is equal to some value at that position
 * @param aBCFlag vector of flag
 * @param iflag_nonzero the value of flag 
 */
template<typename T>
DFM2_INLINE void setRHS_Zero(
    std::vector<T> &vec_b,
    const std::vector<int> &aBCFlag,
    int iflag_nonzero);

std::complex<double> MultSumX(
    const std::complex<double> *va,
    const std::complex<double> *vb,
    unsigned int n);

void XPlusAYBZ(
    std::vector<double> &X,
    size_t nDoF,
    const std::vector<int> &aBCFlag,
    double alpha,
    const std::vector<double> &Y,
    double beta,
    const std::vector<double> &Z);

void XPlusAYBZCW(
    std::vector<double> &X,
    const size_t nDoF,
    const std::vector<int> &aBCFlag,
    double alpha,
    const std::vector<double> &Y,
    double beta,
    const std::vector<double> &Z,
    double gamma,
    const std::vector<double> &W);

DFM2_INLINE void NormalizeX(
    double *p0,
    size_t n);

DFM2_INLINE void OrthogonalizeToUnitVectorX(
    double *p1,
    const double *p0,
    size_t n);

// set boundary condition

void MatVec(
    double* y,
    const double* A, unsigned int ncol, unsigned int nrow,
    const double* x);

/**
 * @param y vector size of nrow
 * @param A a row-major matrix with size [ncol, nrow]
 * @param x vector size of ncol
 */
void MatTVec(
    double* y,
    const double* A, unsigned int ncol, unsigned int nrow,
    const double* x);




// ----------------------------------------------------------------------

/**
 * @brief solve complex linear system using conjugate gradient method
 */
template<typename REAL, class MAT>
std::vector<REAL>
Solve_CG_Complex(
    std::vector<std::complex<REAL>> &r_vec,
    std::vector<std::complex<REAL>> &u_vec,
    REAL conv_ratio_tol,
    unsigned int max_iteration,
    const MAT& mat)
{
  using COMPLEX = std::complex<REAL>;
  assert(!mat.valDia.empty());
  assert(mat.nblk_col == mat.nblk_row);
  assert(mat.nrowdim == mat.len_row);
  const unsigned int ndof = mat.nblk_col * mat.nrowdim;
  assert(r_vec.size() == ndof);
  std::vector<double> aConv;
  u_vec.assign(ndof, 0.0);   // {x} = 0
  double sqnorm_res = Dot(r_vec, r_vec).real();
  if (sqnorm_res < 1.0e-30) { return aConv; }
  const double inv_sqnorm_res_ini = 1.0 / sqnorm_res;
  COMPLEX Ap_vec(ndof);
  COMPLEX p_vec = r_vec;// {p} = {r} (Set Initial Serch Direction)
  for (unsigned int iitr = 0; iitr < max_iteration; iitr++) {
    double alpha;
    {  // alpha = (r,r) / (p,Ap)
      mat.MatVec(Ap_vec.data(),
                 1.0, p_vec.data(), 0.0);
      COMPLEX C_pAp = Dot(p_vec, Ap_vec);
      assert(fabs(C_pAp.imag()) < 1.0e-3);
      const double pAp = C_pAp.real();
      alpha = sqnorm_res / pAp;
    }
    AXPY(COMPLEX(+alpha), p_vec, u_vec);    // {x} = +alpha*{ p} + {x} (updatex)
    AXPY(COMPLEX(-alpha), Ap_vec, r_vec);  // {r} = -alpha*{Ap} + {r}
    double sqnorm_res_new = Dot(r_vec, r_vec).real();
    const double conv_ratio = sqrt(sqnorm_res * inv_sqnorm_res_ini);
    aConv.push_back(conv_ratio);
    if (conv_ratio < conv_ratio_tol) { return aConv; }
    { // update p
      const double beta = sqnorm_res_new / sqnorm_res;  // beta = (r1,r1) / (r0,r0)
      sqnorm_res = sqnorm_res_new;
      for (unsigned int i = 0; i < ndof; i++) { p_vec[i] = r_vec[i] + beta * p_vec[i]; } // {p} = {r} + beta*{p}
    }
  }
  return aConv;
}

template<typename REAL, class MAT>
std::vector<REAL> Solve_BiCGStab(
    std::vector<REAL> &r_vec,
    std::vector<REAL> &x_vec,
    REAL conv_ratio_tol,
    unsigned int max_niter,
    const MAT& mat)
{
  assert(!mat.val_dia_.empty());
  assert(mat.nrowblk_ == mat.ncolblk_);
  assert(mat.nrowdim_ == mat.ncoldim_);
  const unsigned int ndof = mat.nrowblk_ * mat.nrowdim_;
  assert(r_vec.size() == ndof);
  
  std::vector<double> aConv;
  double sq_inv_norm_res_ini;
  {
    const double sq_norm_res_ini = Dot(r_vec, r_vec);
    if (sq_norm_res_ini < 1.0e-30) { return aConv; }
    sq_inv_norm_res_ini = 1.0 / sq_norm_res_ini;
  }
  
  std::vector<double> s_vec(ndof);
  std::vector<double> As_vec(ndof);
  std::vector<double> p_vec(ndof);
  std::vector<double> Ap_vec(ndof);
  std::vector<double> r2_vec(ndof);
  
  x_vec.assign(ndof, 0.0);
  
  r2_vec = r_vec;   // {r2} = {r}
  p_vec = r_vec;    // {p} = {r}
  double r_r2 = Dot(r_vec, r2_vec);   // calc ({r},{r2})
  
  for (unsigned int iitr = 0; iitr < max_niter; iitr++) {
    mat.MatVec(Ap_vec.data(),
               1.0, p_vec.data(), 0.0); // calc {Ap} = [A]*{p}
    double alpha;
    { // alhpa = ({r},{r2}) / ({Ap},{r2})
      const double denominator = Dot(Ap_vec, r2_vec);
      alpha = r_r2 / denominator;
    }
    // {s} = {r} - alpha*{Ap}
    s_vec = r_vec;
    AXPY(-alpha, Ap_vec, s_vec);
    // calc {As} = [A]*{s}
    mat.MatVec(As_vec.data(),
               1.0, s_vec.data(), 0.0);
    // calc omega
    double omega;
    { // omega = ({As},{s}) / ({As},{As})
      const double denominator = Dot(As_vec, As_vec);
      const double numerator = Dot(As_vec, s_vec);
      omega = numerator / denominator;
    }
    // ix += alpha*{p} + omega*{s} (update solution)
    AXPY(alpha, p_vec, x_vec);
    AXPY(omega, s_vec, x_vec);
    // {r} = {s} - omega*{As} (update residual)
    r_vec = s_vec;
    AXPY(-omega, As_vec, r_vec);
    {
      const double sq_norm_res = Dot(r_vec, r_vec);
      const double conv_ratio = sqrt(sq_norm_res * sq_inv_norm_res_ini);
      aConv.push_back(conv_ratio);
      if (conv_ratio < conv_ratio_tol) { return aConv; }
    }
    { // compute beta
      const double tmp1 = Dot(r_vec, r2_vec);
      const double beta = (tmp1 * alpha) / (r_r2 * omega);     // beta = ({r},{r2})^new/({r},{r2})^old * alpha / omega
      r_r2 = tmp1;
      // {p} = {r} + beta*({p}-omega*[A]*{p})  (update p_vector)
      for (unsigned int i = 0; i < ndof; ++i) { p_vec[i] *= beta; }
      AXPY(1.0, r_vec, p_vec);
      AXPY(-beta * omega, Ap_vec, p_vec);
    }
  }
  return aConv;
}

template<typename REAL, class MAT>
std::vector<REAL>
Solve_BiCGSTAB_Complex(
    std::vector<std::complex<REAL>> &r_vec,
    std::vector<std::complex<REAL>> &x_vec,
    REAL conv_ratio_tol,
    unsigned int max_niter,
    const MAT& mat)
{
  using COMPLEX = std::complex<REAL>;
  //
  assert(!mat.valDia.empty());
  assert(mat.nblk_col == mat.nblk_row);
  assert(mat.nrowdim == mat.len_row);
  const unsigned int ndof = mat.nblk_col * mat.nrowdim;
  assert(r_vec.size() == ndof);
  
  std::vector<double> aConv;
  double sq_inv_norm_res_ini;
  {
    const double sq_norm_res_ini = Dot(r_vec, r_vec).real();
    if (sq_norm_res_ini < 1.0e-30) { return aConv; }
    sq_inv_norm_res_ini = 1.0 / sq_norm_res_ini;
  }
  
  std::vector<COMPLEX> s_vec(ndof);
  std::vector<COMPLEX> As_vec(ndof);
  std::vector<COMPLEX> p_vec(ndof);
  std::vector<COMPLEX> Ap_vec(ndof);
  
  x_vec.assign(ndof, 0.0);
  
  const std::vector<COMPLEX> r0_vec = r_vec;   // {r2} = {r}
  p_vec = r_vec;    // {p} = {r}
  COMPLEX r_r0 = Dot(r_vec, r0_vec);   // calc ({r},{r2})
  
  for (unsigned int iitr = 0; iitr < max_niter; iitr++) {
    mat.MatVec(Ap_vec.data(),
               1.0, p_vec.data(), 0.0); // calc {Ap} = [A]*{p}
    const COMPLEX alpha = r_r0 / Dot(Ap_vec, r0_vec); // alhpa = ({r},{r2}) / ({Ap},{r2})
                                                            // {s} = {r} - alpha*{Ap}
    s_vec = r_vec;
    AXPY(-alpha, Ap_vec, s_vec);
    // calc {As} = [A]*{s}
    mat.MatVec(As_vec.data(),
               1.0, s_vec.data(), 0.0);
    // calc omega
    const COMPLEX omega = Dot(s_vec, As_vec) / Dot(As_vec, As_vec).real();  // omega=({As},{s})/({As},{As})
                                                                                        // ix += alpha*{p} + omega*{s} (update solution)
    AXPY(alpha, p_vec, x_vec);
    AXPY(omega, s_vec, x_vec);
    // {r} = {s} - omega*{As} (update residual)
    r_vec = s_vec;
    AXPY(-omega, As_vec, r_vec);
    {
      const double sq_norm_res = Dot(r_vec, r_vec).real();
      const double conv_ratio = sqrt(sq_norm_res * sq_inv_norm_res_ini);
      aConv.push_back(conv_ratio);
      if (conv_ratio < conv_ratio_tol) { return aConv; }
    }
    { // compute beta
      const COMPLEX tmp1 = Dot(r_vec, r0_vec);
      const COMPLEX beta = (tmp1 * alpha) / (r_r0 * omega); // beta = ({r},{r2})^new/({r},{r2})^old * alpha / omega
      r_r0 = tmp1;
      // {p} = {r} + beta*({p}-omega*[A]*{p})  (update p_vector)
      for (unsigned int i = 0; i < ndof; ++i) { p_vec[i] *= beta; }
      AXPY(COMPLEX(1.0), r_vec, p_vec);
      AXPY(-beta * omega, Ap_vec, p_vec);
    }
  }
  return aConv;
}

template <typename REAL, class MAT, class PREC>
std::vector<double> Solve_PBiCGStab(
 REAL* r_vec,
 REAL* x_vec,
 double conv_ratio_tol,
 unsigned int max_niter,
 const MAT& mat,
 const PREC& ilu)
{
  assert( !mat.val_dia_.empty() );
  assert( mat.nrowblk_ == mat.ncolblk_ );
  assert( mat.nrowdim_ == mat.ncoldim_ );
  const unsigned int ndof = mat.nrowblk_*mat.nrowdim_;
  std::vector<double> aResHistry;
  
  // {u} = 0
  for(unsigned int i=0;i<ndof;++i){ x_vec[i] = 0.0; }
  
  double sq_inv_norm_res_ini;
  {
    const double sq_norm_res_ini = DotX(r_vec,r_vec,ndof);
    if( sq_norm_res_ini < 1.0e-60 ){
      aResHistry.push_back( sqrt( sq_norm_res_ini ) );
      return aResHistry;
    }
    sq_inv_norm_res_ini = 1.0 / sq_norm_res_ini;
  }
  
  //    std::cout << "SqIniRes : " << ls.DOT(ir,ir) << std::endl;
  
  std::vector<double> s_vec(ndof);
  std::vector<double> Ms_vec(ndof);
  std::vector<double> AMs_vec(ndof);
  std::vector<double> Mp_vec(ndof);
  std::vector<double> AMp_vec(ndof);
  
  const std::vector<double> r0_vec(r_vec,r_vec+ndof);   // {r2} = {r}
  std::vector<double> p_vec(r_vec,r_vec+ndof);  // {p} = {r}
  
  for(unsigned int iitr=1;iitr<max_niter;iitr++){
    // {Mp_vec} = [M^-1]*{p}
    Mp_vec = p_vec;
    ilu.SolvePrecond(Mp_vec.data());
    // calc (r,r0*)
    const double r_r2 = DotX(r_vec,r0_vec.data(),ndof);
    // calc {AMp_vec} = [A]*{Mp_vec}
    mat.MatVec(AMp_vec.data(),
               1.0, Mp_vec.data(), 0.0);
    // calc alpha
    const double alpha = r_r2 / Dot(AMp_vec,r0_vec);
    // calc s_vector
    s_vec.assign(r_vec,r_vec+ndof);
    AXPY(-alpha,AMp_vec,s_vec);
    // {Ms_vec} = [M^-1]*{s}
    Ms_vec = s_vec;
    ilu.SolvePrecond(Ms_vec.data());
    // calc {AMs_vec} = [A]*{Ms_vec}
    mat.MatVec(AMs_vec.data(),
               1.0,Ms_vec.data(),0.0);
    double omega;
    {  // calc omega
      const double denominator = Dot(AMs_vec,AMs_vec);
      const double numerator = Dot(s_vec,AMs_vec);
      omega = numerator / denominator;
    }
    AXPY(alpha,Mp_vec.data(),x_vec,ndof);
    AXPY(omega,Ms_vec.data(),x_vec,ndof);
    for(unsigned int i=0;i<ndof;++i){ r_vec[i] = s_vec[i]; } // update residual
    AXPY(-omega,AMs_vec.data(),r_vec,ndof);
    {
      const double sq_norm_res = DotX(r_vec,r_vec,ndof);
      const double conv_ratio = sqrt(sq_norm_res * sq_inv_norm_res_ini);
      aResHistry.push_back( conv_ratio );
      if( conv_ratio < conv_ratio_tol ){ return aResHistry; }
    }
    double beta;
    {  // calc beta
      const double tmp1 = DotX(r_vec,r0_vec.data(),ndof);
      beta = tmp1 * alpha / (r_r2*omega);
    }
    // update p_vector
    for(unsigned int i=0;i<ndof;++i){ p_vec[i] *= beta; }
    AXPY(1.0,r_vec,p_vec.data(),ndof);
    AXPY(-beta*omega,AMp_vec,p_vec);
  }
  
  return aResHistry;
}

template <typename REAL, class MAT, class PREC>
std::vector<double> Solve_PBiCGStab_Complex(
    std::complex<REAL>* r_vec,
    std::complex<REAL>* x_vec,
    double conv_ratio_tol,
    unsigned int max_niter,
    const MAT& mat,
    const PREC& ilu)
{
  using COMPLEX = std::complex<REAL>;
  
  assert( !mat.valDia.empty() );
  assert( mat.nblk_col == mat.nblk_row );
  assert( mat.nrowdim == mat.len_row );
  const unsigned int ndof = mat.nblk_col*mat.nrowdim;
  std::vector<double> aResHistry;
  
  for(unsigned int i=0;i<ndof;++i){ x_vec[i] = COMPLEX(0.0,0.0); }   // {u} = 0
  
  double sq_inv_norm_res_ini;
  {
    const double sq_norm_res_ini = Dot(r_vec,r_vec,ndof).real();
    if( sq_norm_res_ini < 1.0e-60 ){
      aResHistry.push_back( sqrt( sq_norm_res_ini ) );
      return aResHistry;
    }
    sq_inv_norm_res_ini = 1.0 / sq_norm_res_ini;
  }
  
  std::vector<COMPLEX> s_vec(ndof);
  std::vector<COMPLEX> Ms_vec(ndof);
  std::vector<COMPLEX> AMs_vec(ndof);
  std::vector<COMPLEX> Mp_vec(ndof);
  std::vector<COMPLEX> AMp_vec(ndof);
  
  const std::vector<COMPLEX> r0_vec(r_vec,r_vec+ndof);   // {r2} = {r}
  std::vector<COMPLEX> p_vec(r_vec,r_vec+ndof);  // {p} = {r}
  
  // calc (r,r0*)
  COMPLEX r_r0 = Dot(r_vec,r0_vec.data(),ndof);
  
  for(unsigned int itr=0;itr<max_niter;itr++){
    // {Mp_vec} = [M^-1]*{p}
    Mp_vec.assign(p_vec.begin(),p_vec.end());
    ilu.SolvePrecond(Mp_vec);
    // calc {AMp_vec} = [A]*{Mp_vec}
    mat.MatVec(AMp_vec.data(),
               COMPLEX(1,0), Mp_vec.data(), COMPLEX(0,0));
    // calc alpha
    const COMPLEX alpha = r_r0 / Dot(AMp_vec,r0_vec);
    // calc s_vector
    s_vec.assign(r_vec,r_vec+ndof);
    AXPY(-alpha,AMp_vec,s_vec);
    // {Ms_vec} = [M^-1]*{s}
    Ms_vec.assign(s_vec.begin(),s_vec.end());
    ilu.SolvePrecond(Ms_vec);
    // calc {AMs_vec} = [A]*{Ms_vec}
    mat.MatVec(AMs_vec.data(),
               COMPLEX(1,0),Ms_vec.data(), COMPLEX(0,0));
    const COMPLEX omega = Dot(s_vec,AMs_vec) / Dot(AMs_vec,AMs_vec).real();
    for(unsigned int i=0;i<ndof;++i){ x_vec[i] = x_vec[i]+alpha*Mp_vec[i]+omega*Ms_vec[i]; }
    for(unsigned int i=0;i<ndof;++i){ r_vec[i] = s_vec[i]-omega*AMs_vec[i]; }
    {
      const double sq_norm_res = Dot(r_vec,r_vec,ndof).real();
      const double conv_ratio = sqrt(sq_norm_res * sq_inv_norm_res_ini);
      aResHistry.push_back( conv_ratio );
      if( conv_ratio < conv_ratio_tol ){ return aResHistry; }
    }
    COMPLEX beta;
    {  // calc beta
      const COMPLEX tmp1 = Dot(r_vec,r0_vec.data(),ndof);
      beta = (tmp1*alpha)/(r_r0*omega);
      r_r0 = tmp1;
    }
    // update p_vector
    for(unsigned int i=0;i<ndof;++i){ p_vec[i] = r_vec[i]+beta*(p_vec[i]-omega*AMp_vec[i]); }
  }
  
  return aResHistry;
}

template <typename REAL, class MAT, class PREC>
std::vector<double> Solve_PCG_Complex(
    std::complex<REAL> *r_vec,
    std::complex<REAL> *x_vec,
    double conv_ratio_tol,
    unsigned int max_nitr,
    const MAT &mat,
    const PREC &ilu)
{
  using COMPLEX = std::complex<REAL>;

  const unsigned int ndof = mat.nblk_col * mat.nrowdim;
  std::vector<double> aResHistry;
  
  for (unsigned int i = 0; i < ndof; i++) { x_vec[i] = COMPLEX(0.0, 0.0); }    // {x} = 0
  
  double inv_sqnorm_res0;
  {
    const double sqnorm_res0 = Dot(r_vec, r_vec, ndof).real();
    aResHistry.push_back(sqnorm_res0);
    if (sqnorm_res0 < 1.0e-30) { return aResHistry; }
    inv_sqnorm_res0 = 1.0 / sqnorm_res0;
  }
  
  // {Pr} = [P]{r}
  std::vector<COMPLEX> Pr_vec(r_vec, r_vec + ndof);
  ilu.SolvePrecond(Pr_vec);
  // {p} = {Pr}
  std::vector<COMPLEX> p_vec = Pr_vec;
  // rPr = ({r},{Pr})
  COMPLEX rPr = Dot(r_vec, Pr_vec.data(), ndof);
  for (unsigned int iitr = 0; iitr < max_nitr; iitr++) {
    {
      std::vector<COMPLEX> &Ap_vec = Pr_vec;
      // {Ap} = [A]{p}
      mat.MatVec(Ap_vec.data(),
                 1.0, p_vec.data(), 0.0);
      // alpha = ({r},{Pr})/({p},{Ap})
      const double pAp = Dot(p_vec, Ap_vec).real();
      COMPLEX alpha = rPr / pAp;
      AXPY(-alpha, Ap_vec.data(), r_vec, ndof);       // {r} = -alpha*{Ap} + {r}
      AXPY(+alpha, p_vec.data(), x_vec, ndof);       // {x} = +alpha*{p } + {x}
    }
    {  // Converge Judgement
      double sqnorm_res = Dot(r_vec, r_vec, ndof).real();
      double conv_ratio = sqrt(sqnorm_res * inv_sqnorm_res0);
      aResHistry.push_back(conv_ratio);
      if (conv_ratio < conv_ratio_tol) { return aResHistry; }
    }
    {  // calc beta
       // {Pr} = [P]{r}
      for (unsigned int i = 0; i < ndof; i++) { Pr_vec[i] = r_vec[i]; }
      ilu.SolvePrecond(Pr_vec);
      // rPr1 = ({r},{Pr})
      const COMPLEX rPr1 = Dot(r_vec, Pr_vec.data(), ndof);
      // beta = rPr1/rPr
      COMPLEX beta = rPr1 / rPr;
      rPr = rPr1;
      // {p} = {Pr} + beta*{p}
      for (unsigned int i = 0; i < ndof; i++) { p_vec[i] = Pr_vec[i] + beta * p_vec[i]; }
    }
  }
  {
    // Converge Judgement
    double sq_norm_res = Dot(r_vec, r_vec, ndof).real();
    aResHistry.push_back(sqrt(sq_norm_res));
  }
  return aResHistry;
}

template <typename REAL, class MAT, class PREC>
std::vector<double> Solve_PCOCG(
    std::complex<REAL>* r_vec,
    std::complex<REAL>* x_vec,
    double conv_ratio_tol,
    unsigned int max_niter,
    const MAT& mat,
    const PREC& ilu)
{
  using COMPLEX = std::complex<REAL>;

  assert( !mat.val_dia_.empty() );
  assert( mat.nrowblk_ == mat.ncolblk_ );
  assert( mat.nrowdim_ == mat.ncoldim_ );
  const unsigned int ndof = mat.nrowblk_*mat.nrowdim_;
  std::vector<double> aResHistry;
  
  for(unsigned int i=0;i<ndof;++i){ x_vec[i] = COMPLEX(0.0,0.0); }   // {u} = 0
  
  double sq_inv_norm_res_ini;
  {
    const double sq_norm_res_ini = DotX(r_vec,r_vec,ndof).real();
    if( sq_norm_res_ini < 1.0e-60 ){
      aResHistry.push_back( sqrt( sq_norm_res_ini ) );
      return aResHistry;
    }
    sq_inv_norm_res_ini = 1.0 / sq_norm_res_ini;
  }
  
  std::vector<COMPLEX> Ap_vec(ndof);
  std::vector<COMPLEX> w_vec(r_vec,r_vec+ndof);
  ilu.SolvePrecond(w_vec.data());
  
  std::vector<COMPLEX> p_vec = w_vec;  // {p} = {w}
  COMPLEX r_w = MultSumX(r_vec,w_vec.data(),ndof);
  
  for(unsigned int itr=0;itr<max_niter;itr++){
    mat.MatVec(Ap_vec.data(),
               COMPLEX(1,0), p_vec.data(), COMPLEX(0,0));
    const COMPLEX alpha = r_w / MultSumX(p_vec.data(),Ap_vec.data(),ndof);
    AXPY(+alpha,p_vec.data(), x_vec,ndof);
    AXPY(-alpha,Ap_vec.data(), r_vec,ndof);
    {
      const double sq_norm_res = DotX(r_vec,r_vec,ndof).real();
      const double conv_ratio = sqrt(sq_norm_res * sq_inv_norm_res_ini);
      aResHistry.push_back( conv_ratio );
      if( conv_ratio < conv_ratio_tol ){ return aResHistry; }
    }
    w_vec.assign(r_vec,r_vec+ndof);
    ilu.SolvePrecond(w_vec.data());
    COMPLEX beta;
    {  // calc beta
      const COMPLEX tmp1 = MultSumX(r_vec,w_vec.data(),ndof);
      beta = tmp1/r_w;
      r_w = tmp1;
    }
    for(unsigned int i=0;i<ndof;++i){ p_vec[i] = w_vec[i] + beta*p_vec[i]; }
  }
  return aResHistry;
}

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/vecxitrsol.cpp"
#endif
  
#endif // MATDIA_CRS_H
