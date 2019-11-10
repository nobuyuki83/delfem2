/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_MATRIX_SPARSE_H
#define DFM2_MATRIX_SPARSE_H

#include <vector>
#include <cassert>
#include <complex>

namespace delfem2 {

template <typename T>
class CMatrixSparse
{
public:
  CMatrixSparse(): nblk_col(0), nblk_row(0), len_col(0), len_row(0) {}
  virtual ~CMatrixSparse(){
    colInd.clear();
    rowPtr.clear();
    valCrs.clear();
    valDia.clear();
  }  
  void Initialize(int nblk, int len, bool is_dia){
    this->nblk_col = nblk;
    this->len_col = len;
    this->nblk_row = nblk;
    this->len_row = len;
    colInd.assign(nblk+1,0);
    rowPtr.clear();
    valCrs.clear();
    if( is_dia ){ valDia.assign(nblk*len*len,0.0); }
    else{         valDia.clear(); }
  }
  void operator = (const CMatrixSparse& m){
    this->nblk_col = m.nblk_col;
    this->len_col  = m.len_col;
    this->nblk_row = m.nblk_row;
    this->len_row  = m.len_row;
    colInd = m.colInd;
    rowPtr = m.rowPtr;
    valCrs = m.valCrs;
    valDia = m.valDia; // copy value
  }
  
  void SetPattern(const int* colind, unsigned int ncolind,
                  const int* rowptr, unsigned int nrowptr){
    assert( rowPtr.empty() );
    assert( ncolind == nblk_col+1 );
    for(unsigned int iblk=0;iblk<nblk_col+1;iblk++){ colInd[iblk] = colind[iblk]; }
    const unsigned int ncrs = colind[nblk_col];
    assert( ncrs == nrowptr );
    rowPtr.resize(ncrs);
    for(unsigned int icrs=0;icrs<ncrs;icrs++){ rowPtr[icrs] = rowptr[icrs]; }
    valCrs.resize(ncrs*len_col*len_row);
  }
  bool SetZero(){
    if( valDia.size() != 0 ){
      assert( len_col == len_row );
      assert( nblk_col == nblk_row );
      const unsigned int n = valDia.size();
      assert( n == len_col*len_col*nblk_col );
      for(unsigned int i=0;i<n;++i){ valDia[i] = 0; }
    }
    {
      const unsigned int n = valCrs.size();
      assert( n == len_col*len_row*rowPtr.size() );
      for(unsigned int i=0;i<n;i++){ valCrs[i] = 0.0; }
    }
    return true;
  }
	bool Mearge(unsigned int nblkel_col, const unsigned int* blkel_col,
              unsigned int nblkel_row, const unsigned int* blkel_row,
              unsigned int blksize, const T* emat,
              std::vector<int>& m_marge_tmp_buffer);
  /**
   * @brief Matrix vector product as: {y} = alpha * [A]{x} + beta * {y}
   */
	void MatVec(T alpha, const std::vector<T>& x, T beta,
              std::vector<T>& y) const;
  /**
   * @brief if BCFlag is 0 for a dof, set all the off-diagonal componenet to zero and set diagonal to one.
   */
  void SetBoundaryCondition(const int* pBCFlag, unsigned int nP, unsigned int ndimVal);
  void AddDia(T eps){
    assert( this->nblk_row == this->nblk_col );
    assert( this->len_row == this->len_col );
    const int blksize = len_col*len_row;
    const int nlen = this->len_col;
    if( valDia.empty() ){ return; }
    for(unsigned int ino=0;ino<nblk_col;++ino){
      for(int ilen=0;ilen<nlen;++ilen){
        valDia[ino*blksize+ilen*nlen+ilen] += eps;
      }
    }
  }
  /**
   * @brief add vector to diagonal component
   * @param[in] lm a lumped mass vector with size of nblk
   * @param[in] scale scaling factor for the lumped mass (typically 1/dt^2).
   * @details the matrix need to be square matrix
   */
  void AddDia_LumpedMass(const T* lm, double scale){
    assert( this->nblk_row == this->nblk_col );
    assert( this->len_row == this->len_col );
    const int blksize = len_col*len_row;
    const int nlen = this->len_col;
    if( valDia.empty() ){ return; }
    for(unsigned int iblk=0;iblk<nblk_col;++iblk){
      for(int ilen=0;ilen<nlen;++ilen){
        valDia[iblk*blksize+ilen*nlen+ilen] += lm[iblk];
      }
    }
  }
public:
	unsigned int nblk_col;
  unsigned int nblk_row;
  unsigned int len_col;
  unsigned int len_row;
  std::vector<unsigned int> colInd;
  std::vector<unsigned int> rowPtr;
  std::vector<T> valCrs;
  std::vector<T> valDia;
};
  
  
// Calc Matrix Vector Product
// {y} = alpha*[A]{x} + beta*{y}
template <>
void CMatrixSparse<double>::MatVec(double alpha,
                                   const std::vector<double>& x,
                                   double beta,
                                   std::vector<double>& y) const;
  
} // end namespace delfem2


// ---------------------------------------------------------------
// implementation of the template functions from here

template<typename T>
void delfem2::CMatrixSparse<T>::MatVec
(T alpha,
 const std::vector<T>& x,
 T beta,
 std::vector<T>& y) const
{
  assert(y.size()==len_col*nblk_col);
  assert(x.size()==len_row*nblk_row);
  const int blksize = len_col*len_row;
  const T* vcrs  = valCrs.data();
  const T* vdia = valDia.data();
  const unsigned int* colind = colInd.data();
  const unsigned int* rowptr = rowPtr.data();
  ////////////////
  for(unsigned int iblk=0;iblk<nblk_col;iblk++){
    for(unsigned int idof=0;idof<len_col;idof++){ y[iblk*len_col+idof] *= beta; }
    const unsigned int colind0 = colind[iblk];
    const unsigned int colind1 = colind[iblk+1];
    for(unsigned int icrs=colind0;icrs<colind1;icrs++){
      assert( icrs < rowPtr.size() );
      const unsigned int jblk0 = rowptr[icrs];
      assert( jblk0 < nblk_row );
      for(unsigned int idof=0;idof<len_col;idof++){
        for(unsigned int jdof=0;jdof<len_row;jdof++){
          y[iblk*len_col+idof] += alpha * vcrs[icrs*blksize+idof*len_col+jdof] * x[jblk0*len_row+jdof];
        }
      }
    }
    for(unsigned int idof=0;idof<len_col;idof++){
      for(unsigned int jdof=0;jdof<len_row;jdof++){
        y[iblk*len_col+idof] += alpha * vdia[iblk*blksize+idof*len_col+jdof] * x[iblk*len_row+jdof];
      }
    }
  }
}

template<typename T>
void delfem2::CMatrixSparse<T>::SetBoundaryCondition
(const int* bc_flag, unsigned int np, unsigned int ndimval)
{
  assert( !this->valDia.empty() );
  assert( this->nblk_row == this->nblk_col );
  assert( this->len_row == this->len_col );
  assert( np == nblk_col );
  assert( ndimval == len_col );
  ////
  const int blksize = len_col*len_row;
  for(unsigned int iblk=0;iblk<nblk_col;iblk++){ // set diagonal
    for(unsigned int ilen=0;ilen<len_col;ilen++){
      if( bc_flag[iblk*len_col+ilen] == 0 ) continue;
      for(unsigned int jlen=0;jlen<len_row;jlen++){
        valDia[iblk*blksize+ilen*len_col+jlen] = 0.0;
        valDia[iblk*blksize+jlen*len_col+ilen] = 0.0;
      }
      valDia[iblk*blksize+ilen*len_col+ilen] = 1.0;
    }
  }
  /////
  for(unsigned int iblk=0;iblk<nblk_col;iblk++){ // set row
    for(unsigned int icrs=colInd[iblk];icrs<colInd[iblk+1];icrs++){
      for(unsigned int ilen=0;ilen<len_col;ilen++){
        if( bc_flag[iblk*len_col+ilen] == 0 ) continue;
        for(unsigned int jlen=0;jlen<len_row;jlen++){
          valCrs[icrs*blksize+ilen*len_col+jlen] = 0.0;
        }
      }
    }
  }
  // -------------
  for(unsigned int icrs=0;icrs<rowPtr.size();icrs++){ // set column
    const int jblk1 = rowPtr[icrs];
    for(unsigned int jlen=0;jlen<len_row;jlen++){
      if( bc_flag[jblk1*len_row+jlen] == 0 ) continue;
      for(unsigned int ilen=0;ilen<len_col;ilen++){
        valCrs[icrs*blksize+ilen*len_col+jlen] = 0.0;
      }
    }
  }
}

template <typename T>
bool delfem2::CMatrixSparse<T>::Mearge
(unsigned int nblkel_col, const unsigned int* blkel_col,
 unsigned int nblkel_row, const unsigned int* blkel_row,
 unsigned int blksize, const T* emat,
 std::vector<int>& marge_buffer)
{
  assert( !valCrs.empty() );
  assert( !valDia.empty() );
  assert( blksize == len_col*len_row );
  marge_buffer.resize(nblk_row);
  const unsigned int* colind = colInd.data();
  const unsigned int* rowptr = rowPtr.data();
  T* vcrs = valCrs.data();
  T* vdia = valDia.data();
  for(unsigned int iblkel=0;iblkel<nblkel_col;iblkel++){
    const unsigned int iblk1 = blkel_col[iblkel];
    assert( iblk1 < nblk_col );
    for(unsigned int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
      assert( jpsup < rowPtr.size() );
      const int jblk1 = rowptr[jpsup];
      marge_buffer[jblk1] = jpsup;
    }
    for(unsigned int jblkel=0;jblkel<nblkel_row;jblkel++){
      const unsigned int jblk1 = blkel_row[jblkel];
      assert( jblk1 < nblk_row );
      if( iblk1 == jblk1 ){  // Marge Diagonal
        const T* pval_in = &emat[(iblkel*nblkel_row+iblkel)*blksize];
        T* pval_out = &vdia[iblk1*blksize];
        for(unsigned int i=0;i<blksize;i++){ pval_out[i] += pval_in[i]; }
      }
      else{  // Marge Non-Diagonal
        if( marge_buffer[jblk1] == -1 ) continue;
        assert( marge_buffer[jblk1] >= 0 && marge_buffer[jblk1] < (int)rowPtr.size() );
        const int jpsup1 = marge_buffer[jblk1];
        assert( rowPtr[jpsup1] == jblk1 );
        const T* pval_in = &emat[(iblkel*nblkel_row+jblkel)*blksize];
        T* pval_out = &vcrs[jpsup1*blksize];
        for(unsigned int i=0;i<blksize;i++){ pval_out[i] += pval_in[i]; }
      }
    }
    for(unsigned int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
      assert( jpsup < rowPtr.size() );
      const int jblk1 = rowptr[jpsup];
      marge_buffer[jblk1] = -1;
    }
  }
  return true;
}


// --------------------------------------------------------------


double CheckSymmetry(const delfem2::CMatrixSparse<double>& mat);
void SetMasterSlave(delfem2::CMatrixSparse<double>& mat, const int* aMSFlag);
void MatSparse_ScaleBlk_LeftRight(delfem2::CMatrixSparse<double>& mat,
                                     const double* scale);
void MatSparse_ScaleBlkLen_LeftRight(delfem2::CMatrixSparse<double>& mat,
                                     const double* scale);




template <typename T>
void XPlusAY(std::vector<T>& X,
             const int nDoF,
             const std::vector<int>& aBCFlag,
             T alpha,
             const std::vector<T>& Y);

template <typename T>
T Dot(const std::vector<T>& r_vec,
      const std::vector<T>& u_vec);

template <typename T>
void AXPY(T a,
          const std::vector<T>& x,
          std::vector<T>& y);

template <typename T>
void AXPY(T a,
          const T* x,
          T* y,
          int n);

template <typename T>
void setRHS_Zero(std::vector<T>& vec_b,
                 const std::vector<int>& aBCFlag,
                 int iflag_nonzero);

template <typename T>
T DotX(const T* r_vec,
       const T* u_vec,
       int ndof);

std::complex<double> MultSumX(const std::complex<double>* va,
                              const std::complex<double>* vb,
                              int n);

void XPlusAYBZ(std::vector<double>& X,
               const int nDoF,
               const std::vector<int>& aBCFlag,
               double alpha,
               const std::vector<double>& Y,
               double beta,
               const std::vector<double>& Z);

void XPlusAYBZCW(std::vector<double>& X,
                 const int nDoF,
                 const std::vector<int>& aBCFlag,
                 double alpha,
                 const std::vector<double>& Y,
                 double beta,
                 const std::vector<double>& Z,
                 double gamma,
                 const std::vector<double>& W);

void ScaleX(double* p0, int n, double s);
void NormalizeX(double* p0, int n);
void OrthogonalizeToUnitVectorX(double* p1,
                                const double* p0, int n);

// set boundary condition

void setRHS_MasterSlave(double* vec_b,
                        int nDoF,
                        const int* aMSFlag);

template <typename T>
std::vector<double>
Solve_CG(std::vector<T>& r_vec,
         std::vector<T>& u_vec,
         double conv_ratio,
         unsigned int iteration,
         const delfem2::CMatrixSparse<T>& mat);

template <typename T>
std::vector<double>
Solve_BiCGSTAB(std::vector<T>& r_vec,
               std::vector<T>& x_vec,
               double conv_ratio,
               unsigned int num_iter,
               const delfem2::CMatrixSparse<T>& mat);

#endif // MATDIA_CRS_H
