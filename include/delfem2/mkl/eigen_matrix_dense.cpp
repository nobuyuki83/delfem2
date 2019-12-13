/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdio.h>
#include <fstream>

#include "delfem2/eigen_densemat.h"

#include "Eigen/Core"
#include "Eigen/SVD"

#ifdef USE_MKL
#include <mkl_cblas.h>
#include <mkl_vml.h>
#include <mkl_lapack.h>
#endif

static double squaredNorm(const std::vector<double>& v)
{
  const int n = (int)v.size();
  double s = 0;
  for (int i = 0; i<n; ++i){ s += v[i]*v[i]; }
  return s;
}

// vai*conj(vbi)
template <typename T>
static double dot
(const std::vector<T>& va, const std::vector<T>& vb)
{
  const int n = (int)va.size();
  double s = 0;
  for (int i = 0; i<n; ++i){ s += va[i]*vb[i]; }
  return s;
}

static void normalize(std::vector<double>& v)
{
  const int n = (int)v.size();
  double sqlen = squaredNorm(v);
  double leninv = 1.0/sqrt(sqlen);
  for (int i = 0; i<n; ++i){ v[i] *= leninv; }
}

// {y} = [A]{x}
template <typename T>
static void matVec
(std::vector<T>& y,
 const std::vector<T>& A,
 const std::vector<T>& x)
{
  const int n = (int)x.size();
  assert(A.size()==n*n);
  y.resize(n);
  for (int i = 0; i<n; ++i){
    double s = 0;
    for (int j = 0; j<n; ++j){ s += A[i*n+j]*x[j]; }
    y[i] = s;
  }
}

// ------------------------------------------------------------
// Solve Matrix with BiCGSTAB Methods
// ------------------------------------------------------------
bool Solve_BiCGSTAB
(double& conv_ratio, int& iteration,
 std::vector<double>& u_vec,
 const std::vector<double>& A,
 const std::vector<double>& y_vec)
{
  
  //	std::cout.precision(18);
  
  const double conv_ratio_tol = conv_ratio;
  const int mx_iter = iteration;
  
  const int n = (int)y_vec.size();
  
  u_vec.assign(n, 0);
  
  assert(A.size()==n*n);
  
  std::vector<double>  r_vec = y_vec;
  std::vector<double>  s_vec(n);
  std::vector<double> As_vec(n);
  std::vector<double>  p_vec(n);
  std::vector<double> Ap_vec(n);
  
  std::vector<double>  r0(n);
  
  double sq_inv_norm_res;
  {
    double dtmp1 = squaredNorm(r_vec);
    //    std::cout << "Initial Residual: " << sqrt(dtmp1) << std::endl;
    if (dtmp1 < 1.0e-30){
      conv_ratio = 0.0;
      iteration = 0;
      return true;
    }
    sq_inv_norm_res = 1.0/dtmp1;
  }
  
  r0 = r_vec;
  //  for (int i = 0; i<n; ++i){ r0[i] = std::conj(r_vec[i]); }
  
  // {p} = {r}
  p_vec = r_vec;
  
  // calc (r,r0*)
  double r_r0conj = dot(r_vec, r0);
  
  iteration = mx_iter;
  for (int iitr = 1; iitr<mx_iter; iitr++){
    
    // calc {Ap} = [A]*{p}
    matVec(Ap_vec, A, p_vec);
    
    // calc alpha
    double alpha;
    {
      const double den = dot(Ap_vec, r0);
      alpha = r_r0conj/den;
    }
    
    // calc s_vector
    for (int i = 0; i<n; ++i){ s_vec[i] = r_vec[i]-alpha*Ap_vec[i]; }
    
    // calc {As} = [A]*{s}
    matVec(As_vec, A, s_vec);
    
    // calc omega
    double omega;
    {
      const double den = squaredNorm(As_vec);
      const double num = dot(As_vec, s_vec);
      omega = num/den;
    }
    
    // update solution
    for (int i = 0; i<n; ++i){ u_vec[i] += alpha*p_vec[i]+omega*s_vec[i]; }
    
    // update residual
    for (int i = 0; i<n; ++i){ r_vec[i] = s_vec[i]-omega*As_vec[i]; }
    
    {
      const double sq_norm_res = squaredNorm(r_vec);
      const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res;
      //      std::cout << iitr << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
      if (sq_conv_ratio < conv_ratio_tol * conv_ratio_tol){
        conv_ratio = sqrt(sq_norm_res * sq_inv_norm_res);
        iteration = iitr;
        return true;
      }
    }
    
    // calc beta
    double beta;
    {
      const double tmp1 = dot(r_vec, r0);
      beta = (tmp1/r_r0conj) * (alpha/omega);
      r_r0conj = tmp1;
    }
    
    // update p_vector
    for (int i = 0; i<n; ++i){
      p_vec[i] = beta*p_vec[i]+r_vec[i]-(beta*omega)*Ap_vec[i];
    }
  }
  return true;
}

template <typename T>
void Solve_CG
(double& conv_ratio,
 int& iteration,  
 std::vector<T>& u_vec,
 const std::vector<T>& A,
 const std::vector<T>& y_vec)
{
  const int n = (int)y_vec.size();
  
  assert(A.size()==n*n);
  
  const double conv_ratio_tol = conv_ratio;
  const int mx_iter = iteration;
  
  std::vector<T>  r_vec = y_vec;
  
  // {x} = 0
  u_vec.assign(n,0.0);
  
  double sqnorm_res = dot(r_vec,r_vec);
  if( sqnorm_res < 1.0e-30 ){
    conv_ratio = 0.0;
    iteration = 0;
    return;
  }
  double inv_sqnorm_res_ini = 1.0 / sqnorm_res;
  
  std::vector<T> Ap_vec(n);
  
  // Set Initial Serch Direction
  // {p} = {r}
  std::vector<T>  p_vec = r_vec;
  
  iteration = mx_iter;
  for(int iitr=1;iitr<mx_iter;iitr++){
    
    double alpha;
    {	// alpha = (r,r) / (p,Ap)
      matVec(Ap_vec, A, p_vec);
      const double pAp = dot(p_vec,Ap_vec);
      alpha = sqnorm_res / pAp;
    }
    
    // update x
    // {x} = +alpha*{ p} + {x}
    for (int i = 0; i<n; ++i){ u_vec[i] += alpha*p_vec[i]; }
    
    // {r} = -alpha*{Ap} + {r}
    for (int i = 0; i<n; ++i){ r_vec[i] -= alpha*Ap_vec[i]; }
    
    double sqnorm_res_new = dot(r_vec,r_vec);
    // Converge Judgement
    
    if( sqnorm_res_new * inv_sqnorm_res_ini < conv_ratio_tol*conv_ratio_tol ){
      conv_ratio = sqrt( sqnorm_res * inv_sqnorm_res_ini );
      iteration = iitr;
      return;
    }
    
    // beta = (r1,r1) / (r0,r0)
    const double beta = sqnorm_res_new / sqnorm_res;
    sqnorm_res = sqnorm_res_new;
    
    // {p} = {r} + beta*{p}
    for(int i=0;i<n;i++){ p_vec[i] = r_vec[i] + beta*p_vec[i]; }
  }
  
  return;
}

////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SolveMatrix
(std::vector<T>& vecR,
 int N,
 std::vector<T>& matA)
{
  std::vector<T> vecTmp = vecR;
  double conv_ratio = 1.0e-5;
  int iteration = 3000;
  //      Solve_BiCGSTAB(conv_ratio,iteration, vecU,matA,vecR);
  Solve_CG(conv_ratio,iteration, vecR,matA,vecTmp);
}  

#ifdef USE_MKL
template <>
void SolveMatrix
(std::vector<double>& vecR,
 int N,
 std::vector<double>& matA)
{
  const int nrhs = 1;
  std::vector<int> iPiv(N,-1);
  int info;
  dgesv(&N,&nrhs, &(matA[0]),&N, &(iPiv[0]), &(vecR[0]),&N, &info );
}
#endif

#ifdef USE_MKL
template <>
void SolveMatrix
(std::vector<float>& vecR,
 int N,
 std::vector<float>& matA)
{
  const int nrhs = 1;
  std::vector<int> iPiv(N,-1);
  int info;
  sgesv(&N,&nrhs, &(matA[0]),&N, &(iPiv[0]), &(vecR[0]),&N, &info );
}
#endif

////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void Subtract
(std::vector<T>& D,
 int n,
 const std::vector<T>& A,
 const typename std::vector<T>::const_iterator B )
{
  for(int k=0;k<n;++k){ D[k] = A[k]-B[k]; }
}

#ifdef USE_MKL
template <>
void Subtract
(std::vector<double>& D,
 int n,
 const std::vector<double>& A,
 const std::vector<double>::const_iterator B )
{
  vdSub(n, &(A[0]), &(B[0]), &(D[0]));
}
#endif

#ifdef USE_MKL
template <>
void Subtract
(std::vector<float>& D,
 int n,
 const std::vector<float>& A,
 const std::vector<float>::const_iterator B )
{
  vsSub(n, &(A[0]), &(B[0]), &(D[0]));
}
#endif

////////////////////////////////////////////////////////////////////////////////////

template <typename T>
double Norm
(int n,
 const std::vector<T>& D)
{
  double W = 0;
  for(int k=0;k<n;++k){
    W += D[k]*D[k];
  }
  return W;
}

#ifdef USE_MKL
template <>
double Norm
(int n,
 const std::vector<double>& D)
{
  return cblas_ddot (n, &(D[0]), 1, &(D[0]), 1);
}
#endif

#ifdef USE_MKL
template <>
double Norm
(int n,
 const std::vector<float>& D)
{
  return cblas_sdot (n, &(D[0]), 1, &(D[0]), 1);
}
#endif

////////////////////////////////////////////////////////////////////////////////////


double myRand(){
  return 2.0*((double)rand()/(RAND_MAX+1.0))-1.0;
  //  return (double)rand()/(RAND_MAX+1.0);
}


void CompSpeed(){
#ifdef USE_MKL
  int N = 1000;
  int M = 1000;
  int nitr = 1000;
  {
    std::vector<float> W(N*M); for(int i=0;i<N*M;++i){ W[ i] = myRand(); }
    std::vector<float> S1(N);  for(int i=0;i<N;  ++i){ S1[i] = myRand(); }
    std::vector<float> aX(M);  for(int i=0;i<M;  ++i){ aX[i] = myRand(); }
    {
      time_t t0 = clock();
      for(int i=0;i<nitr;++i){
        cblas_sgemv(CblasRowMajor,CblasNoTrans, N, M,
                    1.0, &(W[0]), M,
                    &(aX[0]), 1,
                    1.0, &(S1[0]), 1);
      }
      time_t t1 = clock();
      printf("mkl_float: it takes %.2f seconds\n",(double)(t1-t0)/CLOCKS_PER_SEC);
    }

    {
      time_t t0 = clock();
      for(int i=0;i<nitr;++i){
        cblas_sger(CblasRowMajor, N, M,
                   0.1, &(S1[0]), 1,
                   &(aX[0]), 1,
                   &(W[0]), M );
      }
      time_t t1 = clock();
      printf("mkl_float: it takes %.2f seconds\n",(double)(t1-t0)/CLOCKS_PER_SEC);
    }

  }
  {
    std::vector<double> W(N*M); for(int i=0;i<N*M;++i){ W[ i] = myRand(); }
    std::vector<double> S1(N);  for(int i=0;i<N;  ++i){ S1[i] = myRand(); }
    std::vector<double> aX(M);  for(int i=0;i<M;  ++i){ aX[i] = myRand(); }
    {
      time_t t0 = clock();
      for(int i=0;i<nitr;++i){
        cblas_dgemv(CblasRowMajor,CblasNoTrans, N, M,
                    1.0, &(W[0]), M,
                    &(aX[0]), 1,
                    1.0, &(S1[0]), 1);
      }
      time_t t1 = clock();
      printf("double it takes %.2f seconds\n",(double)(t1-t0)/CLOCKS_PER_SEC);
    }

    {
      time_t t0 = clock();
      for(int i=0;i<nitr;++i){
        cblas_dger(CblasRowMajor, N, M,
                   0.1, &(S1[0]), 1,
                   &(aX[0]), 1,
                   &(W[0]), M );
      }
      time_t t1 = clock();
      printf("double: it takes %.2f seconds\n",(double)(t1-t0)/CLOCKS_PER_SEC);
    }
  }
#endif
}

template <typename T>
void shuffleInOut
(std::vector<T>& X1,
 std::vector<T>& R1,
 int ntraining,
 std::vector<int>& ShuffleIndex,
 int ndimIn, int ndimOut,
 const std::vector<T>& aX,
 const std::vector<T>& aR)
{
  for(int it=0;it<ntraining;++it){
    int jt = ShuffleIndex[it];
    for(int i=0;i<ndimIn; i++){ X1[jt*ndimIn +i] = aX[it*ndimIn +i]; }
    for(int j=0;j<ndimOut;j++){ R1[jt*ndimOut+j] = aR[it*ndimOut+j]; }
  }
}


#ifdef USE_MKL
template <>
void shuffleInOut
(std::vector<double>& X1,
 std::vector<double>& R1,
 int ntraining,
 std::vector<int>& ShuffleIndex,
 int ndimIn, int ndimOut,
 const std::vector<double>& aX,
 const std::vector<double>& aR)
{
  for(int it=0;it<ntraining;++it){
    int jt = ShuffleIndex[it];
    cblas_dcopy(ndimIn,  &(aX[it*ndimIn ]), 1, &(X1[jt*ndimIn ]), 1);
    cblas_dcopy(ndimOut, &(aR[it*ndimOut]), 1, &(R1[jt*ndimOut]), 1);
  }
}
#endif


#ifdef USE_MKL
template <>
void shuffleInOut
(std::vector<float>& X1,
 std::vector<float>& R1,
 int ntraining,
 std::vector<int>& ShuffleIndex,
 int ndimIn, int ndimOut,
 const std::vector<float>& aX,
 const std::vector<float>& aR)
{
  for(int it=0;it<ntraining;++it){
    int jt = ShuffleIndex[it];
    cblas_scopy(ndimIn,  &(aX[it*ndimIn ]), 1, &(X1[jt*ndimIn ]), 1);
    cblas_scopy(ndimOut, &(aR[it*ndimOut]), 1, &(R1[jt*ndimOut]), 1);
  }
}
#endif


void Shuffle(std::vector<int>& X)
{
  const int N = (int)X.size();
  for(int i=0;i<N;i++){
    int j = rand()%N;
    int k = X[i];
    X[i] = X[j];
    X[j] = k;
  }
}


////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void AddDiagonal
(std::vector<T>& A,
 int N, double alpha, double beta)
{
  double sum = 0;
  for(int i=0;i<N;++i){
    sum += fabs(A[i*N+i]);
  }
  sum /= N;
  for(int i=0;i<N;++i){
    A[i*N+i] = A[i*N+i]*alpha;
    A[i*N+i] += sum*beta;
  }
}

template <typename T>
void MergeLMA
(std::vector<T>& A,
 std::vector<T>& u_vec,
 int N, int M,
 const std::vector<T>& J,
 const std::vector<T>& FR)
{
  for(int k=0;k<M;++k){
    for(int i=0;i<N;++i){
      for(int j=0;j<N;++j){
        A[i*N+j] += J[k*N+i]*J[k*N+j];
      }
    }
  }
  for(int i=0;i<N;++i){
    for(int k=0;k<M;++k){
      u_vec[i] -= FR[k]*J[k*N+i];
    }
  }
}

#ifdef USE_MKL
template <>
void MergeLMA
(std::vector<double>& A,
 std::vector<double>& u_vec,
 int N, int M,
 const std::vector<double>& J,
 const std::vector<double>& FR)
{
  {
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                N, N, M, 1.0,
                &(J[0]), N,
                &(J[0]), N,
                1.0, &(A[0]), N);
  }
  for(int i=0;i<N;++i){
    for(int k=0;k<M;++k){
      u_vec[i] -= FR[k]*J[k*N+i];
    }
  }

}
#endif

#ifdef USE_MKL
template <>
void MergeLMA
(std::vector<float>& A,
 std::vector<float>& u_vec,
 int N, int M,
 const std::vector<float>& J,
 const std::vector<float>& FR)
{
  {
    cblas_sgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                N, N, M, 1.0,
                &(J[0]), N,
                &(J[0]), N,
                1.0, &(A[0]), N);
  }
  for(int i=0;i<N;++i){
    for(int k=0;k<M;++k){
      u_vec[i] -= FR[k]*J[k*N+i];
    }
  }
}
#endif


////////////////////////////////////////////////////////////////////

template <typename T>
void WXPB
(std::vector<T>& S1,
 const int nt,
 const typename std::vector<T>::const_iterator aX,
 const int ndimIn, const int ndimOut,
 const std::vector<T>& W,
 const std::vector<T>& b)
{
  assert( S1.size() >= nt*ndimOut );
  for(int it=0;it<nt;++it){
    for(int j=0;j<ndimOut;++j){
      double uj = b[j];
      for(int i=0;i<ndimIn;++i){
        uj += W[j*ndimIn+i]*aX[it*ndimIn+i];
      }
      S1[it*ndimOut+j] = uj;
    }
  }
}

#ifdef USE_MKL
// double
template <>
void WXPB
(std::vector<double>& S1,
 const int nt,
 const typename std::vector<double>::const_iterator aX,
 const int ndimIn, const int ndimOut,
 const std::vector<double>& W,
 const std::vector<double>& b)
{
  assert( S1.size() >= nt*ndimOut );
  for(int it=0;it<nt;++it){
    cblas_dcopy(ndimOut, &(b[0]), 1, &(S1[it*ndimOut]), 1);
  }
  if( nt == 1 ){
    cblas_dgemv(CblasRowMajor,CblasNoTrans, ndimOut, ndimIn,
                1.0, &(W[0]), ndimIn,
                &(aX[0]), 1,
                1.0, &(S1[0]), 1);
  }
  else{
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                nt, ndimOut, ndimIn, 1.0,
                &(aX[0]), ndimIn,
                &(W[0]), ndimIn,
                1.0, &(S1[0]), ndimOut);
  }
}

//float
template <>
void WXPB
(std::vector<float>& S1,
 const int nt,
 const typename std::vector<float>::const_iterator aX,
 const int ndimIn, const int ndimOut,
 const std::vector<float>& W,
 const std::vector<float>& b)
{
  assert( S1.size() >= nt*ndimOut );
  for(int it=0;it<nt;++it){
    cblas_scopy(ndimOut, &(b[0]), 1, &(S1[it*ndimOut]), 1);
  }
  if( nt == 1 ){
    cblas_sgemv(CblasRowMajor,CblasNoTrans, ndimOut, ndimIn,
                1.0, &(W[0]), ndimIn,
                &(aX[0]), 1,
                1.0, &(S1[0]), 1);
  }
  else{
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                nt, ndimOut, ndimIn, 1.0,
                &(aX[0]), ndimIn,
                &(W[0]), ndimIn,
                1.0, &(S1[0]), ndimOut);
  }
}
#endif


///////////////////////////////////////////////////////////////////////////////////

template <typename T>
void WTX
(typename std::vector<T>::iterator D0,
 int nt, const typename std::vector<T>::const_iterator D1,
 int ndimIn, int ndimOut,
 const std::vector<T>& W)
{
  assert( W.size() == ndimIn*ndimOut );
  for(int it=0;it<nt;++it){
    for(int i=0;i<ndimIn;++i){
      double dj = 0;
      for(int j=0;j<ndimOut;++j){
        dj += D1[it*ndimOut+j]*W[j*ndimIn+i];
      }
      D0[it*ndimIn+i] = dj;
    }
  }
}

#ifdef USE_MKL
// double
template <>
void WTX
(typename std::vector<double>::iterator D0,
 int nt, const typename std::vector<double>::const_iterator D1,
 int ndimIn, int ndimOut,
 const std::vector<double>& W)
{
  assert( W.size() == ndimIn*ndimOut );
  if( nt == 1 ){
    cblas_dgemv(CblasRowMajor,CblasTrans, ndimOut, ndimIn,
                1.0, &(W[0]), ndimIn,
                &(D1[0]), 1,
                0.0, &(D0[0]), 1);
  }
  else{
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                nt, ndimIn, ndimOut, 1.0,
                &(D1[0]), ndimOut,
                &(W[0]), ndimIn,
                0.0, &(D0[0]), ndimIn);
  }
}

// float
template <>
void WTX
(typename std::vector<float>::iterator D0,
 int nt, const typename std::vector<float>::const_iterator D1,
 int ndimIn, int ndimOut,
 const std::vector<float>& W)
{
  assert( W.size() == ndimIn*ndimOut );
  if( nt == 1 ){
    cblas_sgemv(CblasRowMajor,CblasTrans, ndimOut, ndimIn,
                1.0, &(W[0]), ndimIn,
                &(D1[0]), 1,
                0.0, &(D0[0]), 1);
  }
  else{
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                nt, ndimIn, ndimOut, 1.0,
                &(D1[0]), ndimOut,
                &(W[0]), ndimIn,
                0.0, &(D0[0]), ndimIn);
  }
}
#endif

////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void UpdateWB
(std::vector<T>& W,
 std::vector<T>& b,
 double alpha,
 int ndimIn, int ndimOut, int nt,
 const std::vector<T>& D1,
 const typename std::vector<T>::const_iterator& S0)
{
  for(int it=0;it<nt;++it){
    for(int j=0;j<ndimOut;++j){
      for(int i=0;i<ndimIn;++i){
        W[j*ndimIn+i] += alpha*D1[j]*S0[i];
      }
      b[j] += alpha*D1[j];
    }
  }
}

#ifdef USE_MKL
// double
template <>
void UpdateWB
(std::vector<double>& W,
 std::vector<double>& b,
 double alpha,
 int ndimIn, int ndimOut, int nt,
 const std::vector<double>& D1,
 const typename std::vector<double>::const_iterator& S0)
{
  if( nt == 1 ){
    cblas_daxpy(ndimOut,
                +alpha, &(D1[0]), 1,
                &(b[0]), 1);
    cblas_dger(CblasRowMajor, ndimOut, ndimIn,
               +alpha, &(D1[0]), 1,
               &(S0[0]), 1,
               &(W[0]), ndimIn );
  }
  else{
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                ndimOut, ndimIn, nt, +alpha,
                &(D1[0]), ndimOut,
                &(S0[0]), ndimIn,
                1.0, &(W[0]), ndimIn);
    for(int it=0;it<nt;++it){
      for(int j=0;j<ndimOut;++j){
        b[j] += alpha*D1[it*ndimOut+j];
      }
    }
  }
}

// float
template <>
void UpdateWB
(std::vector<float>& W,
 std::vector<float>& b,
 double alpha,
 int ndimIn, int ndimOut, int nt,
 const std::vector<float>& D1,
 const typename std::vector<float>::const_iterator& S0)
{
  if( nt == 1 ){
    cblas_saxpy(ndimOut,
                +alpha, &(D1[0]), 1,
                &(b[0]), 1);
    cblas_sger(CblasRowMajor, ndimOut, ndimIn,
               +alpha, &(D1[0]), 1,
               &(S0[0]), 1,
               &(W[0]), ndimIn );
  }
  else{
    cblas_sgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                ndimOut, ndimIn, nt, +alpha,
                &(D1[0]), ndimOut,
                &(S0[0]), ndimIn,
                1.0, &(W[0]), ndimIn);
    for(int it=0;it<nt;++it){
      for(int j=0;j<ndimOut;++j){
        b[j] += alpha*D1[it*ndimOut+j];
      }
    }
  }
}
#endif




