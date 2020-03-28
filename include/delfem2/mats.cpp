/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <cmath>
#include <vector>
#include <complex>
#include "delfem2/mats.h"

typedef std::complex<double> COMPLEX;
namespace dfm2 = delfem2;


static double MatNorm_Assym(
                     const double* V0,
                     unsigned int n0,
                     unsigned int m0,
                     const double* V1)
{
  double s = 0.0;
  for(unsigned int i=0;i<n0;++i){
    for(unsigned int j=0;j<m0;++j){
      double v0 = V0[i*m0+j];
      double v1 = V1[j*n0+i];
      s += (v0-v1)*(v0-v1);
    }
  }
  return s;
}

static double MatNorm(
               const double* V,
               unsigned int n,
               unsigned int m)
{
  double s = 0.0;
  for(unsigned int i=0;i<n;++i){
    for(unsigned int j=0;j<m;++j){
      double v = V[i*m+j];
      s += v*v;
    }
  }
  return s;
}

static double MatNorm_Assym(
                     const double* V,
                     unsigned int n)
{
  double s = 0.0;
  for(unsigned int i=0;i<n;++i){
    for(unsigned int j=0;j<n;++j){
      double v0 = V[i*n+j];
      double v1 = V[j*n+i];
      s += (v0-v1)*(v0-v1);
    }
  }
  return s;
}


// -------------------------------------------------------

// Calc Matrix Vector Product
// {y} = alpha*[A]{x} + beta*{y}
template <typename T>
void dfm2::CMatrixSparse<T>::MatVec
(T* y,
 T alpha,
 const T* x,
 T beta) const
{
  const unsigned int ndofcol = len_col*nblk_col;
  for(unsigned int i=0;i<ndofcol;++i){ y[i] *= beta; }
  const unsigned int blksize = len_col*len_row;
  // --------
	if( len_col == 1 && len_row == 1 ){
		const T* vcrs  = valCrs.data();
		const T* vdia = valDia.data();
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		//
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			T& vy = y[iblk];
			const unsigned int colind0 = colind[iblk];
			const unsigned int colind1 = colind[iblk+1];
			for(unsigned int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < rowPtr.size() );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
				vy += alpha * vcrs[icrs] * x[jblk0];
			}
			vy += alpha * vdia[iblk] * x[iblk];
		}
	}
	else if( len_col == 2 && len_row == 2 ){
		const T* vcrs  = valCrs.data();
		const T* vdia = valDia.data();
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		//
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < rowPtr.size() );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
				y[iblk*2+0] += alpha * ( vcrs[icrs*4  ]*x[jblk0*2+0] + vcrs[icrs*4+1]*x[jblk0*2+1] );
				y[iblk*2+1] += alpha * ( vcrs[icrs*4+2]*x[jblk0*2+0] + vcrs[icrs*4+3]*x[jblk0*2+1] );
			}
			y[iblk*2+0] += alpha * ( vdia[iblk*4+0]*x[iblk*2+0] + vdia[iblk*4+1]*x[iblk*2+1] );
			y[iblk*2+1] += alpha * ( vdia[iblk*4+2]*x[iblk*2+0] + vdia[iblk*4+3]*x[iblk*2+1] );
		}
	}
	else if( len_col == 3 && len_row == 3 ){
		const T* vcrs  = valCrs.data();
		const T* vdia = valDia.data();
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		//
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < rowPtr.size() );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
        const unsigned int i0 = iblk*3;
        const unsigned int j0 = jblk0*3;
        const unsigned int k0 = icrs*9;
				y[i0+0] += alpha*(vcrs[k0+0]*x[j0+0]+vcrs[k0+1]*x[j0+1]+vcrs[k0+2]*x[j0+2]);
				y[i0+1] += alpha*(vcrs[k0+3]*x[j0+0]+vcrs[k0+4]*x[j0+1]+vcrs[k0+5]*x[j0+2]);
				y[i0+2] += alpha*(vcrs[k0+6]*x[j0+0]+vcrs[k0+7]*x[j0+1]+vcrs[k0+8]*x[j0+2]);
			}
      {
        const unsigned int i0 = iblk*3;
        const unsigned int k0 = iblk*9;
        y[i0+0] += alpha*(vdia[k0+0]*x[i0+0]+vdia[k0+1]*x[i0+1]+vdia[k0+2]*x[i0+2]);
        y[i0+1] += alpha*(vdia[k0+3]*x[i0+0]+vdia[k0+4]*x[i0+1]+vdia[k0+5]*x[i0+2]);
        y[i0+2] += alpha*(vdia[k0+6]*x[i0+0]+vdia[k0+7]*x[i0+1]+vdia[k0+8]*x[i0+2]);
      }
		}
  }
	else if( len_col == 4 && len_row == 4 ){
    const T* vcrs  = valCrs.data();
    const T* vdia = valDia.data();
    const unsigned int* colind = colInd.data();
    const unsigned int* rowptr = rowPtr.data();
    //
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
        assert( icrs < rowPtr.size() );
        const unsigned int jblk0 = rowptr[icrs];
        assert( jblk0 < nblk_row );
        const unsigned int i0 = iblk*4;
        const unsigned int j0 = jblk0*4;
        const unsigned int k0 = icrs*16;
        y[i0+0] += alpha*(vcrs[k0+ 0]*x[j0+0]+vcrs[k0+ 1]*x[j0+1]+vcrs[k0+ 2]*x[j0+2]+vcrs[k0+ 3]*x[j0+3]);
        y[i0+1] += alpha*(vcrs[k0+ 4]*x[j0+0]+vcrs[k0+ 5]*x[j0+1]+vcrs[k0+ 6]*x[j0+2]+vcrs[k0+ 7]*x[j0+3]);
        y[i0+2] += alpha*(vcrs[k0+ 8]*x[j0+0]+vcrs[k0+ 9]*x[j0+1]+vcrs[k0+10]*x[j0+2]+vcrs[k0+11]*x[j0+3]);
        y[i0+3] += alpha*(vcrs[k0+12]*x[j0+0]+vcrs[k0+13]*x[j0+1]+vcrs[k0+14]*x[j0+2]+vcrs[k0+15]*x[j0+3]);
      }
      {
        const unsigned int i0 = iblk*4;
        const unsigned int k0 = iblk*16;
        y[i0+0] += alpha*(vdia[k0+ 0]*x[i0+0]+vdia[k0+ 1]*x[i0+1]+vdia[k0+ 2]*x[i0+2]+vdia[k0+ 3]*x[i0+3]);
        y[i0+1] += alpha*(vdia[k0+ 4]*x[i0+0]+vdia[k0+ 5]*x[i0+1]+vdia[k0+ 6]*x[i0+2]+vdia[k0+ 7]*x[i0+3]);
        y[i0+2] += alpha*(vdia[k0+ 8]*x[i0+0]+vdia[k0+ 9]*x[i0+1]+vdia[k0+10]*x[i0+2]+vdia[k0+11]*x[i0+3]);
        y[i0+3] += alpha*(vdia[k0+12]*x[i0+0]+vdia[k0+13]*x[i0+1]+vdia[k0+14]*x[i0+2]+vdia[k0+15]*x[i0+3]);
      }
    }
  }
	else{
		const T* vcrs  = valCrs.data();
		const T* vdia = valDia.data();
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		//
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			const unsigned int colind0 = colind[iblk];
			const unsigned int colind1 = colind[iblk+1];
			for(unsigned int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < rowPtr.size() );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
				for(unsigned int idof=0;idof<len_col;idof++){
          for(unsigned int jdof=0;jdof<len_row;jdof++){
            y[iblk*len_col+idof] += alpha * vcrs[icrs*blksize+idof*len_row+jdof] * x[jblk0*len_row+jdof];
          }
				}
			}
			for(unsigned int idof=0;idof<len_col;idof++){
        for(unsigned int jdof=0;jdof<len_row;jdof++){
          y[iblk*len_col+idof] += alpha * vdia[iblk*blksize+idof*len_row+jdof] * x[iblk*len_row+jdof];
        }
			}
		}
	}
}
template void delfem2::CMatrixSparse<float>::MatVec(float *y, float alpha, const float *x, float beta) const;
template void delfem2::CMatrixSparse<double>::MatVec(double *y, double alpha, const double *x, double beta) const;
template void delfem2::CMatrixSparse<COMPLEX>::MatVec(COMPLEX *y, COMPLEX alpha, const COMPLEX *x, COMPLEX beta) const;

// -------------------------------------------------------

// Calc Matrix Vector Product
// {y} = alpha*[A]^T{x} + beta*{y}
template <typename T>
void dfm2::CMatrixSparse<T>::MatTVec
 (T* y,
  T alpha,
  const T* x,
  T beta) const
{
  const unsigned int ndofrow = len_row*nblk_row;
  for(unsigned int i=0;i<ndofrow;++i){ y[i] *= beta; }
  const unsigned int blksize = len_col*len_row;
  // ---------
  /*
  if( len_col == 1 && len_row == 1 ){
    const T* vcrs  = valCrs.data();
    const T* vdia = valDia.data();
    const unsigned int* colind = colInd.data();
    const unsigned int* rowptr = rowPtr.data();
    //
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
      T& vy = y[iblk];
      vy *= beta;
      const unsigned int colind0 = colind[iblk];
      const unsigned int colind1 = colind[iblk+1];
      for(unsigned int icrs=colind0;icrs<colind1;icrs++){
        assert( icrs < rowPtr.size() );
        const unsigned int jblk0 = rowptr[icrs];
        assert( jblk0 < nblk_row );
        vy += alpha * vcrs[icrs] * x[jblk0];
      }
      vy += alpha * vdia[iblk] * x[iblk];
    }
  }
  else if( len_col == 2 && len_row == 2 ){
    const T* vcrs  = valCrs.data();
    const T* vdia = valDia.data();
    const unsigned int* colind = colInd.data();
    const unsigned int* rowptr = rowPtr.data();
    //
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
      y[iblk*2+0] *= beta;
      y[iblk*2+1] *= beta;
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
        assert( icrs < rowPtr.size() );
        const unsigned int jblk0 = rowptr[icrs];
        assert( jblk0 < nblk_row );
        y[iblk*2+0] += alpha * ( vcrs[icrs*4  ]*x[jblk0*2+0] + vcrs[icrs*4+1]*x[jblk0*2+1] );
        y[iblk*2+1] += alpha * ( vcrs[icrs*4+2]*x[jblk0*2+0] + vcrs[icrs*4+3]*x[jblk0*2+1] );
      }
      y[iblk*2+0] += alpha * ( vdia[iblk*4+0]*x[iblk*2+0] + vdia[iblk*4+1]*x[iblk*2+1] );
      y[iblk*2+1] += alpha * ( vdia[iblk*4+2]*x[iblk*2+0] + vdia[iblk*4+3]*x[iblk*2+1] );
    }
  }
  else if( len_col == 3 && len_row == 3 ){
    const T* vcrs  = valCrs.data();
    const T* vdia = valDia.data();
    const unsigned int* colind = colInd.data();
    const unsigned int* rowptr = rowPtr.data();
    //
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
      y[iblk*3+0] *= beta;
      y[iblk*3+1] *= beta;
      y[iblk*3+2] *= beta;
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
        assert( icrs < rowPtr.size() );
        const unsigned int jblk0 = rowptr[icrs];
        assert( jblk0 < nblk_row );
        const unsigned int i0 = iblk*3;
        const unsigned int j0 = jblk0*3;
        const unsigned int k0 = icrs*9;
        y[i0+0] += alpha*(vcrs[k0+0]*x[j0+0]+vcrs[k0+1]*x[j0+1]+vcrs[k0+2]*x[j0+2]);
        y[i0+1] += alpha*(vcrs[k0+3]*x[j0+0]+vcrs[k0+4]*x[j0+1]+vcrs[k0+5]*x[j0+2]);
        y[i0+2] += alpha*(vcrs[k0+6]*x[j0+0]+vcrs[k0+7]*x[j0+1]+vcrs[k0+8]*x[j0+2]);
      }
      {
        const unsigned int i0 = iblk*3;
        const unsigned int k0 = iblk*9;
        y[i0+0] += alpha*(vdia[k0+0]*x[i0+0]+vdia[k0+1]*x[i0+1]+vdia[k0+2]*x[i0+2]);
        y[i0+1] += alpha*(vdia[k0+3]*x[i0+0]+vdia[k0+4]*x[i0+1]+vdia[k0+5]*x[i0+2]);
        y[i0+2] += alpha*(vdia[k0+6]*x[i0+0]+vdia[k0+7]*x[i0+1]+vdia[k0+8]*x[i0+2]);
      }
    }
  }
  else if( len_col == 4 && len_row == 4 ){
    const T* vcrs  = valCrs.data();
    const T* vdia = valDia.data();
    const unsigned int* colind = colInd.data();
    const unsigned int* rowptr = rowPtr.data();
    //
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
      y[iblk*4+0] *= beta;
      y[iblk*4+1] *= beta;
      y[iblk*4+2] *= beta;
      y[iblk*4+3] *= beta;
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
        assert( icrs < rowPtr.size() );
        const unsigned int jblk0 = rowptr[icrs];
        assert( jblk0 < nblk_row );
        const unsigned int i0 = iblk*4;
        const unsigned int j0 = jblk0*4;
        const unsigned int k0 = icrs*16;
        y[i0+0] += alpha*(vcrs[k0+ 0]*x[j0+0]+vcrs[k0+ 1]*x[j0+1]+vcrs[k0+ 2]*x[j0+2]+vcrs[k0+ 3]*x[j0+3]);
        y[i0+1] += alpha*(vcrs[k0+ 4]*x[j0+0]+vcrs[k0+ 5]*x[j0+1]+vcrs[k0+ 6]*x[j0+2]+vcrs[k0+ 7]*x[j0+3]);
        y[i0+2] += alpha*(vcrs[k0+ 8]*x[j0+0]+vcrs[k0+ 9]*x[j0+1]+vcrs[k0+10]*x[j0+2]+vcrs[k0+11]*x[j0+3]);
        y[i0+3] += alpha*(vcrs[k0+12]*x[j0+0]+vcrs[k0+13]*x[j0+1]+vcrs[k0+14]*x[j0+2]+vcrs[k0+15]*x[j0+3]);
      }
      {
        const unsigned int i0 = iblk*4;
        const unsigned int k0 = iblk*16;
        y[i0+0] += alpha*(vdia[k0+ 0]*x[i0+0]+vdia[k0+ 1]*x[i0+1]+vdia[k0+ 2]*x[i0+2]+vdia[k0+ 3]*x[i0+3]);
        y[i0+1] += alpha*(vdia[k0+ 4]*x[i0+0]+vdia[k0+ 5]*x[i0+1]+vdia[k0+ 6]*x[i0+2]+vdia[k0+ 7]*x[i0+3]);
        y[i0+2] += alpha*(vdia[k0+ 8]*x[i0+0]+vdia[k0+ 9]*x[i0+1]+vdia[k0+10]*x[i0+2]+vdia[k0+11]*x[i0+3]);
        y[i0+3] += alpha*(vdia[k0+12]*x[i0+0]+vdia[k0+13]*x[i0+1]+vdia[k0+14]*x[i0+2]+vdia[k0+15]*x[i0+3]);
      }
    }
  }
  else{
   */
  {
    const T* vcrs  = valCrs.data();
    const T* vdia = valDia.data();
    const unsigned int* colind = colInd.data();
    const unsigned int* rowptr = rowPtr.data();
    //
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
      const unsigned int colind0 = colind[iblk];
      const unsigned int colind1 = colind[iblk+1];
      for(unsigned int icrs=colind0;icrs<colind1;icrs++){
        assert( icrs < rowPtr.size() );
        const unsigned int jblk0 = rowptr[icrs];
        assert( jblk0 < nblk_row );
        for(unsigned int idof=0;idof<len_col;idof++){
          for(unsigned int jdof=0;jdof<len_row;jdof++){
            y[jblk0*len_row+jdof] += alpha * vcrs[icrs*blksize+idof*len_row+jdof] * x[iblk*len_col+idof];
          }
        }
      }
      for(unsigned int jdof=0;jdof<len_row;jdof++){
        for(unsigned int idof=0;idof<len_col;idof++){
          y[iblk*len_row+jdof] += alpha * vdia[iblk*blksize+idof*len_row+jdof] * x[iblk*len_col+idof];
        }
      }
    }
  }
}
template void delfem2::CMatrixSparse<float>::MatTVec(float *y, float alpha, const float *x, float beta) const;
template void delfem2::CMatrixSparse<double>::MatTVec(double *y, double alpha, const double *x, double beta) const;
template void delfem2::CMatrixSparse<COMPLEX>::MatTVec(COMPLEX *y, COMPLEX alpha, const COMPLEX *x, COMPLEX beta) const;


// ----------------------------------

template<typename T>
bool delfem2::CMatrixSparse<T>::Mearge
(unsigned int nblkel_col, const unsigned int *blkel_col,
 unsigned int nblkel_row, const unsigned int *blkel_row,
 unsigned int blksize, const T *emat,
 std::vector<int> &marge_buffer)
{
  assert(!valCrs.empty());
  assert(!valDia.empty());
  assert(blksize == len_col * len_row);
  marge_buffer.resize(nblk_row);
  const unsigned int *colind = colInd.data();
  const unsigned int *rowptr = rowPtr.data();
  T *vcrs = valCrs.data();
  T *vdia = valDia.data();
  for (unsigned int iblkel = 0; iblkel < nblkel_col; iblkel++) {
    const unsigned int iblk1 = blkel_col[iblkel];
    assert(iblk1 < nblk_col);
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < rowPtr.size());
      const int jblk1 = rowptr[jpsup];
      marge_buffer[jblk1] = jpsup;
    }
    for (unsigned int jblkel = 0; jblkel < nblkel_row; jblkel++) {
      const unsigned int jblk1 = blkel_row[jblkel];
      assert(jblk1 < nblk_row);
      if (iblk1 == jblk1) {  // Marge Diagonal
        const T *pval_in = &emat[(iblkel * nblkel_row + iblkel) * blksize];
        T *pval_out = &vdia[iblk1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      }
      else {  // Marge Non-Diagonal
        if (marge_buffer[jblk1] == -1) {
          assert(0);
          return false;
        }
        assert(marge_buffer[jblk1] >= 0 && marge_buffer[jblk1] < (int) rowPtr.size());
        const int jpsup1 = marge_buffer[jblk1];
        assert(rowPtr[jpsup1] == jblk1);
        const T *pval_in = &emat[(iblkel * nblkel_row + jblkel) * blksize];
        T *pval_out = &vcrs[jpsup1 * blksize];
        for (unsigned int i = 0; i < blksize; i++) { pval_out[i] += pval_in[i]; }
      }
    }
    for (unsigned int jpsup = colind[iblk1]; jpsup < colind[iblk1 + 1]; jpsup++) {
      assert(jpsup < rowPtr.size());
      const int jblk1 = rowptr[jpsup];
      marge_buffer[jblk1] = -1;
    }
  }
  return true;
}
template bool delfem2::CMatrixSparse<float>::Mearge(unsigned int nblkel_col, const unsigned int *blkel_col,
                                                    unsigned int nblkel_row, const unsigned int *blkel_row,
                                                    unsigned int blksize, const float *emat,
                                                    std::vector<int> &marge_buffer);
template bool delfem2::CMatrixSparse<double>::Mearge(unsigned int nblkel_col, const unsigned int *blkel_col,
                                                     unsigned int nblkel_row, const unsigned int *blkel_row,
                                                     unsigned int blksize, const double *emat,
                                                     std::vector<int> &marge_buffer);
template bool delfem2::CMatrixSparse<COMPLEX>::Mearge(unsigned int nblkel_col, const unsigned int *blkel_col,
                                                      unsigned int nblkel_row, const unsigned int *blkel_row,
                                                      unsigned int blksize, const COMPLEX *emat,
                                                      std::vector<int> &marge_buffer);

// -----------------------------------------------------------------

template<typename T>
void delfem2::CMatrixSparse<T>::SetFixedBC_Dia(
    const int *bc_flag,
    T val_dia)
{
  assert(!this->valDia.empty());
  assert(this->nblk_row == this->nblk_col);
  assert(this->len_row == this->len_col);
  const int blksize = len_col * len_row;
  for (unsigned int iblk = 0; iblk < nblk_col; iblk++) { // set diagonal
    for (unsigned int ilen = 0; ilen < len_col; ilen++) {
      if (bc_flag[iblk * len_col + ilen] == 0) continue;
      for (unsigned int jlen = 0; jlen < len_row; jlen++) {
        valDia[iblk * blksize + ilen * len_col + jlen] = 0.0;
        valDia[iblk * blksize + jlen * len_col + ilen] = 0.0;
      }
      valDia[iblk * blksize + ilen * len_col + ilen] = val_dia;
    }
  }
}
template void delfem2::CMatrixSparse<float>::SetFixedBC_Dia(const int *bc_flag, float val_dia);
template void delfem2::CMatrixSparse<double>::SetFixedBC_Dia(const int *bc_flag, double val_dia);
template void delfem2::CMatrixSparse<COMPLEX>::SetFixedBC_Dia(const int *bc_flag, COMPLEX val_dia);


template<typename T>
void delfem2::CMatrixSparse<T>::SetFixedBC_Row(const int *bc_flag)
{
  assert(!this->valDia.empty());
  assert(this->nblk_row == this->nblk_col);
  assert(this->len_row == this->len_col);
  const int blksize = len_col * len_row;
  for (unsigned int iblk = 0; iblk < nblk_col; iblk++) { // set row
    for (unsigned int icrs = colInd[iblk]; icrs < colInd[iblk + 1]; icrs++) {
      for (unsigned int ilen = 0; ilen < len_col; ilen++) {
        if (bc_flag[iblk * len_col + ilen] == 0) continue;
        for (unsigned int jlen = 0; jlen < len_row; jlen++) {
          valCrs[icrs * blksize + ilen * len_col + jlen] = 0.0;
        }
      }
    }
  }
}
template void delfem2::CMatrixSparse<float>::SetFixedBC_Row(const int *bc_flag);
template void delfem2::CMatrixSparse<double>::SetFixedBC_Row(const int *bc_flag);
template void delfem2::CMatrixSparse<COMPLEX>::SetFixedBC_Row(const int *bc_flag);

template<typename T>
void delfem2::CMatrixSparse<T>::SetFixedBC_Col(const int *bc_flag)
{
  assert(!this->valDia.empty());
  assert(this->nblk_row == this->nblk_col);
  assert(this->len_row == this->len_col);
  const int blksize = len_col * len_row;
  for (unsigned int icrs = 0; icrs < rowPtr.size(); icrs++) { // set column
    const int jblk1 = rowPtr[icrs];
    for (unsigned int jlen = 0; jlen < len_row; jlen++) {
      if (bc_flag[jblk1 * len_row + jlen] == 0) continue;
      for (unsigned int ilen = 0; ilen < len_col; ilen++) {
        valCrs[icrs * blksize + ilen * len_col + jlen] = 0.0;
      }
    }
  }
}
template void delfem2::CMatrixSparse<float>::SetFixedBC_Col(const int *bc_flag);
template void delfem2::CMatrixSparse<double>::SetFixedBC_Col(const int *bc_flag);
template void delfem2::CMatrixSparse<COMPLEX>::SetFixedBC_Col(const int *bc_flag);

// -----------------------------------------------------------------

void dfm2::SetMasterSlave
(dfm2::CMatrixSparse<double>& mat,
 const int* aMSFlag)
{
  assert( !mat.valDia.empty() );
  assert( mat.nblk_row == mat.nblk_col );
  assert( mat.len_row == mat.len_col );
  const unsigned int len = mat.len_col;
  const unsigned int nblk = mat.nblk_col;
  const unsigned int blksize = len*len;
  const unsigned int ndof = nblk*len;
  /////
  std::vector<int> row2crs(nblk,-1);
  for(unsigned int idof1=0;idof1<ndof;++idof1){ // add row
    int idof0 = aMSFlag[idof1];
    if( idof0 == -1 ) continue;
    unsigned int ino0 = idof0 / len;
    unsigned int ilen0 = idof0 - ino0*len;
    assert( ilen0 < len );
    assert( ino0 < nblk && ilen0 < len );
    unsigned int ino1 = idof1 / len;
    unsigned int ilen1 = idof1 - ino1*len;
    assert( ino1 < nblk && ilen1 < len );
    assert( ilen0 == ilen1 );
    for(unsigned int icrs0=mat.colInd[ino0];icrs0<mat.colInd[ino0+1];++icrs0){
      unsigned int jno0 = mat.rowPtr[icrs0];
      assert( jno0 < nblk );
      row2crs[jno0] = icrs0;
    }
    for(unsigned int icrs1=mat.colInd[ino1];icrs1<mat.colInd[ino1+1];++icrs1){
      unsigned int jno1 = mat.rowPtr[icrs1];
      assert( jno1 < nblk );
      assert( jno1 != ino1 );
      if( jno1 != ino0 ){ // add non-diagonal 1 to non-diagonal 0
        const int icrs0 = row2crs[jno1];
        assert( icrs0 >= 0 && icrs0 < (int)mat.rowPtr.size() );
        for(unsigned int jdim=0;jdim<len;++jdim){
          mat.valCrs[icrs0*blksize+ilen0*len+jdim] += mat.valCrs[icrs1*blksize+ilen1*len+jdim];
        }
      }
      else{ // add non-diagonal 1 to diagonal 0
        for(unsigned int jdim=0;jdim<len;++jdim){
          mat.valDia[ino0*blksize+ilen0*len+jdim] += mat.valCrs[icrs1*blksize+ilen1*len+jdim];
        }
      }
    }
    { // add diagonal 1 to non-diagonal 0
      const int icrs0 = row2crs[ino1];
      assert( icrs0 >= 0 && icrs0 < (int)mat.rowPtr.size() );
      for(unsigned int jdim=0;jdim<len;++jdim){
        mat.valCrs[icrs0*blksize+ilen0*len+jdim] += mat.valDia[ino1*blksize+ilen1*len+jdim];
      }
    }
    for(unsigned int icrs0=mat.colInd[ino0];icrs0<mat.colInd[ino0+1];++icrs0){
      int jno0 = mat.rowPtr[icrs0];
      assert( jno0 >= 0 && jno0 < (int)nblk );
      row2crs[jno0] = -1;
    }
  }
  // ---------------------------------------------
  row2crs.assign(nblk,-1);
  for(unsigned int ino=0;ino<nblk;ino++){
    for(unsigned int icrs=mat.colInd[ino];icrs<mat.colInd[ino+1];++icrs){
      int jno0 = mat.rowPtr[icrs];
      assert( jno0 >= 0 && jno0 < (int)nblk );
      row2crs[jno0] = icrs;
    }
    for(unsigned int jlen1=0;jlen1<len;jlen1++){
      int jdof0 = aMSFlag[ino*len+jlen1];
      if( jdof0 == -1 ) continue;
      int jno0 = (int)(jdof0/len);
      assert( jdof0 - jno0*len == jlen1 );
      const int icrs0 = row2crs[jno0];
      assert( icrs0 >= 0 && icrs0 < (int)mat.rowPtr.size() );
      for(unsigned int ilen=0;ilen<len;ilen++){
        mat.valCrs[icrs0*blksize+ilen*len+jlen1] += mat.valDia[ino*blksize+ilen*len+jlen1];
      }
    }
    for(unsigned int icrs1=mat.colInd[ino];icrs1<mat.colInd[ino+1];icrs1++){
      const unsigned int jno1 = mat.rowPtr[icrs1];
      assert( jno1 < nblk );
      for(unsigned int jlen1=0;jlen1<len;jlen1++){
        if( aMSFlag[jno1*len+jlen1] == -1 ) continue;
        auto jdof0 = (unsigned int)aMSFlag[jno1*len+jlen1];
        unsigned int jno0 = jdof0/len;
        assert( jno0 < nblk );
        assert( jdof0 - jno0*len == jlen1 );
        if( ino == jno0 ){
          for(unsigned int ilen=0;ilen<len;ilen++){
            mat.valDia[jno0*blksize+ilen*len+jlen1] += mat.valCrs[icrs1*blksize+ilen*len+jlen1];
          }
        }
        else{
          const int icrs0 = row2crs[jno0];
          assert( icrs0 >= 0 && icrs0 < (int)mat.rowPtr.size() );
          for(unsigned int ilen=0;ilen<len;ilen++){
            mat.valCrs[icrs0*blksize+ilen*len+jlen1] += mat.valCrs[icrs1*blksize+ilen*len+jlen1];
          }
        }
      }
    }
    for(unsigned int icrs=mat.colInd[ino];icrs<mat.colInd[ino+1];++icrs){
      unsigned int jno0 = mat.rowPtr[icrs];
      assert( jno0 < nblk );
      row2crs[jno0] = -1;
    }
  }
  // --------------------------------------
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ilen=0;ilen<len;ilen++){
      if( aMSFlag[iblk*len+ilen] == -1 ) continue;
      for(unsigned int jlen=0;jlen<len;jlen++){
        mat.valDia[iblk*blksize+ilen*len+jlen] = 0.0;
        mat.valDia[iblk*blksize+jlen*len+ilen] = 0.0;
      }
      mat.valDia[iblk*blksize+ilen*len+ilen] = 1.0;
    }
  }
  // ---------------------------------------------
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int icrs=mat.colInd[iblk];icrs<mat.colInd[iblk+1];icrs++){
      for(unsigned int idim=0;idim<len;idim++){
        if( aMSFlag[iblk*len+idim] == -1 ) continue;
        auto idof0 = (unsigned int)aMSFlag[iblk*len+idim];
        unsigned int jblk = mat.rowPtr[icrs];
        for(unsigned int jdim=0;jdim<len;jdim++){
          unsigned int idof1 = jblk*len+jdim;
          if( idof0 != idof1 ){ mat.valCrs[icrs*blksize+idim*len+jdim] = +0.0; }
          else{                 mat.valCrs[icrs*blksize+idim*len+jdim] = -1.0; }
          mat.valCrs[icrs*blksize+idim*len+jdim] = +0.0;
        }
      }
    }
  }
  // ---------------------------------------------
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int icrs=mat.colInd[iblk];icrs<mat.colInd[iblk+1];icrs++){
      const int jblk1 = mat.rowPtr[icrs];
      for(unsigned int jdim=0;jdim<len;jdim++){
        if( aMSFlag[jblk1*len+jdim] == -1 ) continue;
        auto idof0 = (unsigned int)aMSFlag[jblk1*len+jdim];
        for(unsigned int idim=0;idim<len;idim++){
          unsigned int idof1 = iblk*len+idim;
          if( idof0 != idof1 ){ mat.valCrs[icrs*blksize+idim*len+jdim] = +0.0; }
          else{                 mat.valCrs[icrs*blksize+idim*len+jdim] = -1.0; }
          mat.valCrs[icrs*blksize+idim*len+jdim] = +0.0;
        }
      }
    }
  }
}

void dfm2::MatSparse_ScaleBlk_LeftRight
(dfm2::CMatrixSparse<double>& mat,
 const double* scale)
{
  assert( mat.nblk_row == mat.nblk_col );
  assert( mat.len_row == mat.len_col );
  const unsigned int nblk = mat.nblk_col;
  const unsigned int len = mat.len_col;
  const unsigned int blksize = len*len;
  for(unsigned int ino=0;ino<nblk;++ino){
    for(unsigned int icrs0=mat.colInd[ino];icrs0<mat.colInd[ino+1];++icrs0){
      const int jno = mat.rowPtr[icrs0];
      const double s0 = scale[ino]*scale[jno];
      for(unsigned int i=0;i<blksize;++i){ mat.valCrs[icrs0*blksize+i] *= s0; }
    }
  }
  if( !mat.valDia.empty() ){
    for(unsigned int ino=0;ino<nblk;++ino){
      double s0 = scale[ino]*scale[ino];
      for(unsigned int i=0;i<blksize;++i){ mat.valDia[ino*blksize+i] *= s0; }
    }
  }
}

void dfm2::MatSparse_ScaleBlkLen_LeftRight
(dfm2::CMatrixSparse<double>& mat,
 const double* scale)
{
  assert( mat.nblk_row == mat.nblk_col );
  assert( mat.len_row == mat.len_col );
  const unsigned int nblk = mat.nblk_col;
  const unsigned int len = mat.len_col;
  const unsigned int blksize = len*len;
  for(unsigned int ino=0;ino<nblk;++ino){
    for(unsigned int icrs0=mat.colInd[ino];icrs0<mat.colInd[ino+1];++icrs0){
      const int jno = mat.rowPtr[icrs0];
      for(unsigned int ilen=0;ilen<len;++ilen){
        for(unsigned int jlen=0;jlen<len;++jlen){
          mat.valCrs[icrs0*blksize+ilen*len+jlen] *= scale[ino*len+ilen]*scale[jno*len+jlen];
        }
      }
    }
  }
  if( !mat.valDia.empty() ){
    for(unsigned int ino=0;ino<nblk;++ino){
      for(unsigned int ilen=0;ilen<len;++ilen){
        for(unsigned int jlen=0;jlen<len;++jlen){
          mat.valDia[ino*blksize+ilen*len+jlen] *= scale[ino*len+ilen]*scale[ino*len+jlen];
        }
      }
    }
  }
}

double dfm2::CheckSymmetry(
    const dfm2::CMatrixSparse<double>& mat)
{
  assert( mat.nblk_row == mat.nblk_col );
  assert( mat.len_row == mat.len_col );
  const unsigned int blksize = mat.len_col*mat.len_row;
  const unsigned int nlen = mat.len_col;
  //
  double sum = 0;
  for(unsigned int ino=0;ino<mat.nblk_col;++ino){
    for(unsigned int icrs0=mat.colInd[ino];icrs0<mat.colInd[ino+1];++icrs0){
      int jno = mat.rowPtr[icrs0];
      unsigned int icrs1 = mat.colInd[jno];
      for(;icrs1<mat.colInd[jno+1];++icrs1){
        if( mat.rowPtr[icrs1] == ino ){ break; }
      }
      if( icrs1 == mat.colInd[jno+1] ){ // no counterpart
        sum += MatNorm(mat.valCrs.data()+blksize*icrs0, mat.len_col, mat.len_row);
      }
      else{
        sum += MatNorm_Assym(mat.valCrs.data()+blksize*icrs0, mat.len_col, mat.len_row,
                             mat.valCrs.data()+blksize*icrs1);
      }
    }
    sum += MatNorm_Assym(mat.valDia.data()+blksize*ino,nlen);
  }
  return sum;
}
