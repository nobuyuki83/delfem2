/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cassert>
#include <vector>
#include <complex>

#include "delfem2/lsilu_mats.h"

// ----------------------------------------------------

namespace delfem2 {
namespace ilu {

DFM2_INLINE void CalcMatPr(
    double *out, const double *d, double *tmp,
    const unsigned int ni, const unsigned int nj)
{
  unsigned int i, j, k;
  for (i = 0; i < ni; i++) {
    for (j = 0; j < nj; j++) {
      tmp[i * nj + j] = 0.0;
      for (k = 0; k < ni; k++) {
        tmp[i * nj + j] += d[i * ni + k] * out[k * nj + j];
      }
    }
  }
  for (i = 0; i < ni * nj; i++) {
    out[i] = tmp[i];
  }
}

DFM2_INLINE void CalcSubMatPr(
    double *out, const double *a, const double *b,
    const int ni, const int nk, const int nj)
{
  int i, j, k;
  for (i = 0; i < ni; i++) {
    for (j = 0; j < nj; j++) {
      for (k = 0; k < nk; k++) {
        out[i * nj + j] -= a[i * nk + k] * b[k * nj + j];
      }
    }
  }
}


DFM2_INLINE void CalcInvMat(
    double *a,
    const unsigned int n,
    int &info)
{
  double tmp1;

  info = 0;
  unsigned int i, j, k;
  for (i = 0; i < n; i++) {
    if (fabs(a[i * n + i]) < 1.0e-30) {
      info = 1;
      return;
    }
    if (a[i * n + i] < 0.0) {
      info--;
    }
    tmp1 = 1.0 / a[i * n + i];
    a[i * n + i] = 1.0;
    for (k = 0; k < n; k++) {
      a[i * n + k] *= tmp1;
    }
    for (j = 0; j < n; j++) {
      if (j != i) {
        tmp1 = a[j * n + i];
        a[j * n + i] = 0.0;
        for (k = 0; k < n; k++) {
          a[j * n + k] -= tmp1 * a[i * n + k];
        }
      }
    }
  }
}

// t is a tmporary buffer size of 9
DFM2_INLINE void CalcInvMat3(double a[], double t[])
{
  const double det =
      +a[0] * a[4] * a[8] + a[3] * a[7] * a[2] + a[6] * a[1] * a[5]
      - a[0] * a[7] * a[5] - a[6] * a[4] * a[2] - a[3] * a[1] * a[8];
  const double inv_det = 1.0 / det;

  for (int i = 0; i < 9; i++) { t[i] = a[i]; }

  a[0] = inv_det * (t[4] * t[8] - t[5] * t[7]);
  a[1] = inv_det * (t[2] * t[7] - t[1] * t[8]);
  a[2] = inv_det * (t[1] * t[5] - t[2] * t[4]);

  a[3] = inv_det * (t[5] * t[6] - t[3] * t[8]);
  a[4] = inv_det * (t[0] * t[8] - t[2] * t[6]);
  a[5] = inv_det * (t[2] * t[3] - t[0] * t[5]);

  a[6] = inv_det * (t[3] * t[7] - t[4] * t[6]);
  a[7] = inv_det * (t[1] * t[6] - t[0] * t[7]);
  a[8] = inv_det * (t[0] * t[4] - t[1] * t[3]);
}

class CRowLev {
public:
  CRowLev() : row(0), lev(0) {}

  CRowLev(int row, int lev) : row(row), lev(lev) {}

  bool operator<(const CRowLev &rhs) const {
    if (row != rhs.row) {
      return (row < rhs.row);
    }
    return (lev < rhs.lev);
  }

public:
  int row;
  int lev;
};

class CRowLevNext {
public:
  unsigned int row;
  unsigned int lev;
  int next;
};

} // ilu
} // delfem2


// -----------------------------------------------------------------


template <typename T>
delfem2::CPreconditionerILU<T>::CPreconditionerILU(const CPreconditionerILU<T>& p)
{
//  std::cout << "CPreconditionerILU -- construct copy" << std::endl;
  this->mat = p.mat; // deep copy
  const unsigned int nblk = this->mat.nblk_col;
  this->m_diaInd.resize(nblk);
  for(int iblk=0;iblk<nblk;++iblk){
    this->m_diaInd[iblk] = p.m_diaInd[iblk];
  }
}

// -------------------------------------------------------------------

namespace delfem2{

// numerical factorization
template <>
DFM2_INLINE bool
CPreconditionerILU<double>::DoILUDecomp()
{
  const int nmax_sing = 10;
	int icnt_sing = 0;
  
	const unsigned int len = mat.nrowdim;
  const unsigned int nblk = mat.nrowblk;
//  const int m_ncrs = mat.m_ncrs;
  const unsigned int* colind = mat.colInd.data();
  const unsigned int* rowptr = mat.rowPtr.data();
  double* vcrs = mat.valCrs.data();
  double* vdia = mat.valDia.data();
#ifndef NDEBUG
  const unsigned int m_ncrs = colind[nblk];
#endif
  
  std::vector<int> row2crs(nblk,-1);
  
	if( len == 1 ){
		for(unsigned int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<colind[nblk] );
				const unsigned int jblk0 = rowptr[ijcrs];
        assert( jblk0<nblk );
				row2crs[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U] 
			for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
				const unsigned int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double ikvalue = vcrs[ikcrs];
				for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const unsigned int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs[jblk0];
						if( ijcrs0 == -1 ) continue;
						vcrs[ijcrs0] -= ikvalue*vcrs[kjcrs];
					}
					else{ vdia[iblk] -= ikvalue*vcrs[kjcrs]; }
				}
			}
			double iivalue = vdia[iblk];
			if( fabs(iivalue) > 1.0e-30 ){
				vdia[iblk] = 1.0 / iivalue;
			}
			else{
				std::cout << "frac false" << iblk << std::endl;
				icnt_sing++;
				if( icnt_sing > nmax_sing ){
					return false;
				}
			}
      for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				vcrs[ijcrs] = vcrs[ijcrs] * vdia[iblk];
			}
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const unsigned int jblk0 = rowptr[ijcrs];
        assert( jblk0<nblk );
				row2crs[jblk0] = -1;
			}
		}	// end iblk
	}
	// -----------------------------
	else if( len == 2 ){
		double TmpBlk[4];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const unsigned int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U]
			for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
				const unsigned int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double* vik = &vcrs[ikcrs*4];
				for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const unsigned int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
					double* vkj = &vcrs[kjcrs*4]; assert( vkj != 0 );
					double* vij = 0;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs[jblk0];
						if( ijcrs0 == -1 ) continue;	
						vij = &vcrs[ijcrs0*4];
					}
					else{
            vij = &vdia[iblk*4];
          }
          assert( vij != nullptr );
					vij[0] -= vik[0]*vkj[0]+vik[1]*vkj[2];
					vij[1] -= vik[0]*vkj[1]+vik[1]*vkj[3];
					vij[2] -= vik[2]*vkj[0]+vik[3]*vkj[2];
					vij[3] -= vik[2]*vkj[1]+vik[3]*vkj[3];
				}
			}
			{
				double* vii = &vdia[iblk*4];
				const double det = vii[0]*vii[3]-vii[1]*vii[2];
				if( fabs(det) > 1.0e-30 ){
					const double inv_det = 1.0/det;
					double dtmp1 = vii[0];
					vii[0] =  inv_det*vii[3];
					vii[1] = -inv_det*vii[1];
					vii[2] = -inv_det*vii[2];
					vii[3] =  inv_det*dtmp1;
				}
				else{
					std::cout << "frac false" << iblk << std::endl;
					icnt_sing++;
					if( icnt_sing > nmax_sing ){
						return false;
					}
				}
			}
			// [U] = [1/D][U]
      for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				double* pVal_ij = &vcrs[ijcrs*4];
				const double* vii = &vdia[iblk*4];
				for(int i=0;i<4;i++){ TmpBlk[i] = pVal_ij[i]; }
				pVal_ij[0] = vii[0]*TmpBlk[0] + vii[1]*TmpBlk[2];
				pVal_ij[1] = vii[0]*TmpBlk[1] + vii[1]*TmpBlk[3];
				pVal_ij[2] = vii[2]*TmpBlk[0] + vii[3]*TmpBlk[2];
				pVal_ij[3] = vii[2]*TmpBlk[1] + vii[3]*TmpBlk[3];
			}
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const unsigned int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = -1;
			}
		}	// end iblk
	}
	// -----------------------------------------------------------
	else if( len == 3 ){	// lenBlk >= 3
		double tmpBlk[9];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const unsigned int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U]
			for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
				const unsigned int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double* vik = &vcrs[ikcrs*9];
				for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const unsigned int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
					double* vkj = &vcrs[kjcrs*9]; assert( vkj != 0 );
					double* vij = nullptr;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs[jblk0];
            if( ijcrs0 == -1 ){ continue; }
						vij = &vcrs[ijcrs0*9];
					}
					else{
            vij = &vdia[iblk*9];
          }
					assert( vij != nullptr );
          for(int i=0;i<3;i++){
            vij[i*3+0] -= vik[i*3+0]*vkj[0] + vik[i*3+1]*vkj[3] + vik[i*3+2]*vkj[6];
            vij[i*3+1] -= vik[i*3+0]*vkj[1] + vik[i*3+1]*vkj[4] + vik[i*3+2]*vkj[7];
            vij[i*3+2] -= vik[i*3+0]*vkj[2] + vik[i*3+1]*vkj[5] + vik[i*3+2]*vkj[8];
          }
				}
			}
			{  
				double* vii = &vdia[iblk*9];
				const double det =
        + vii[0]*vii[4]*vii[8] + vii[3]*vii[7]*vii[2] + vii[6]*vii[1]*vii[5]
        - vii[0]*vii[7]*vii[5] - vii[6]*vii[4]*vii[2] - vii[3]*vii[1]*vii[8];
				if( fabs(det) > 1.0e-30 ){
					ilu::CalcInvMat3(vii,tmpBlk);
				}
				else{
					std::cout << "frac false 3 " << iblk << std::endl;
					icnt_sing++;
					if( icnt_sing > nmax_sing ){
            std::cout << "ilu frac false exceeds tolerance" << std::endl;
						return false;
					}
				}
			}
			// [U] = [1/D][U]
      for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				double* vij = &vcrs[ijcrs*9];
				const double* vii = &vdia[iblk*9];
        for(int i=0;i<9;i++){ tmpBlk[i] = vij[i]; }
        for(int i=0;i<3;i++){
          vij[i*3+0] = vii[i*3+0]*tmpBlk[0] + vii[i*3+1]*tmpBlk[3] + vii[i*3+2]*tmpBlk[6];
          vij[i*3+1] = vii[i*3+0]*tmpBlk[1] + vii[i*3+1]*tmpBlk[4] + vii[i*3+2]*tmpBlk[7];
          vij[i*3+2] = vii[i*3+0]*tmpBlk[2] + vii[i*3+1]*tmpBlk[5] + vii[i*3+2]*tmpBlk[8];
        }
			}		
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const unsigned int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = -1;
			}
		}	// end iblk
	}
  // ------------------------------------------------------------------------
	else{	// lenBlk >= 4
    const unsigned int blksize = len*len;
		auto* pTmpBlk = new double [blksize];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const unsigned int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U] 
			for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
				const unsigned int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double* vik = &vcrs[ikcrs*blksize];
				for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const unsigned int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
					double* vkj = &vcrs[kjcrs *blksize]; assert( vkj != 0 );
					double* vij = nullptr;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs[jblk0];
            if( ijcrs0 == -1 ){ continue; }
						vij = &vcrs[ijcrs0*blksize];
					}
					else{
            vij = &vdia[iblk *blksize];
          }
					assert( vij != nullptr );
          ilu::CalcSubMatPr(vij,vik,vkj, len,len,len);
				}
			}
			{
				double* vii = &vdia[iblk*blksize];
				int info = 0;
				ilu::CalcInvMat(vii,len,info);
				if( info==1 ){
					std::cout << "frac false" << iblk << std::endl;
					icnt_sing++;
					if( icnt_sing > nmax_sing ){
						delete[] pTmpBlk;
						return false;
					}
				}
			}
			// [U] = [1/D][U]
      for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				double* vij = &vcrs[ijcrs*blksize];
				const double* pVal_ii = &vdia[iblk *blksize];
        ilu::CalcMatPr(vij,pVal_ii,pTmpBlk,  len,len);
			}
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const unsigned int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = -1;
			}
		}	// end iblk
		delete[] pTmpBlk;
	}  
	return true;
}

// numerical factorization
template <>
DFM2_INLINE bool
CPreconditionerILU<std::complex<double>>::DoILUDecomp()
{
  typedef std::complex<double> COMPLEX;
//  const int nmax_sing = 10;
//  int icnt_sing = 0;
  
  const unsigned int len = mat.nrowdim;
  const unsigned int nblk = mat.nrowblk;
  //  const int m_ncrs = mat.m_ncrs;
  const unsigned int* colind = mat.colInd.data();
  const unsigned int* rowptr = mat.rowPtr.data();
  COMPLEX* vcrs = mat.valCrs.data();
  COMPLEX* vdia = mat.valDia.data();
#ifndef NDEBUG
  const unsigned int m_ncrs = colind[nblk];
#endif
  
  std::vector<int> row2crs(nblk,-1);
  if( len == 1 ){
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<colind[nblk] );
        const unsigned int jblk0 = rowptr[ijcrs];
        assert( jblk0<nblk );
        row2crs[jblk0] = ijcrs;
      }
      // [L] * [D^-1*U]
      for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
        const unsigned int kblk = rowptr[ikcrs]; assert( kblk<nblk );
        const COMPLEX ikvalue = vcrs[ikcrs];
        for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
          const unsigned int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
          if( jblk0 != iblk ){
            const int ijcrs0 = row2crs[jblk0];
            if( ijcrs0 == -1 ) continue;
            vcrs[ijcrs0] -= ikvalue*vcrs[kjcrs];
          }
          else{ vdia[iblk] -= ikvalue*vcrs[kjcrs]; }
        }
      }
      COMPLEX iivalue = vdia[iblk];
      vdia[iblk] = 1.0/iivalue;
      /*
      if( fabs((std::conj(iivalue)*iivalue).real()) > 1.0e-30 ){

      }
      else{
        std::cout << "frac false" << iblk << std::endl;
        icnt_sing++;
        if( icnt_sing > nmax_sing ){
          return false;
        }
      }
       */
      for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
        vcrs[ijcrs] = vcrs[ijcrs] * vdia[iblk];
      }
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
        const unsigned int jblk0 = rowptr[ijcrs];
        assert( jblk0<nblk );
        row2crs[jblk0] = -1;
      }
    }  // end iblk
  }
  else{
    std::cout << "error!-->TOBE IMPLEMENTED" << std::endl;
    abort();
  }
  return true;
}

}


// -----------------------------------------------------

template <typename T>
void delfem2::CPreconditionerILU<T>::ForwardSubstitution
( T* vec ) const
{
  const unsigned int len = mat.nrowdim;
  const unsigned int nblk = mat.nrowblk;
  
  if( len == 1 ){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    const T* vdia = mat.valDia.data();
    // -------------------------
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      T lvec_i = vec[iblk];
      for(unsigned int ijcrs=colind[iblk];ijcrs<m_diaInd[iblk];ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = rowptr[ijcrs];
        assert( jblk0<(int)iblk );
        lvec_i -= vcrs[ijcrs]*vec[jblk0];
      }
      vec[iblk] = vdia[iblk]*lvec_i;
    }
  }
  else if( len == 2 ){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    const T* vdia = mat.valDia.data();
    // ------------------------
    T pTmpVec[2];
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      pTmpVec[0] = vec[iblk*2+0];
      pTmpVec[1] = vec[iblk*2+1];
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = m_diaInd[iblk];
      for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = rowptr[ijcrs];
        assert( jblk0<(int)iblk );
        const T* vij = &vcrs[ijcrs*4];
        const T valj0 = vec[jblk0*2+0];
        const T valj1 = vec[jblk0*2+1];
        pTmpVec[0] -= vij[0]*valj0+vij[1]*valj1;
        pTmpVec[1] -= vij[2]*valj0+vij[3]*valj1;
      }
      const T* vii = &vdia[iblk*4];
      vec[iblk*2+0] = vii[0]*pTmpVec[0]+vii[1]*pTmpVec[1];
      vec[iblk*2+1] = vii[2]*pTmpVec[0]+vii[3]*pTmpVec[1];
    }
  }
  else if( len == 3 ){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    const T* vdia = mat.valDia.data();
    // -------------------------
    T pTmpVec[3];
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      pTmpVec[0] = vec[iblk*3+0];
      pTmpVec[1] = vec[iblk*3+1];
      pTmpVec[2] = vec[iblk*3+2];
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = m_diaInd[iblk];
      for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = rowptr[ijcrs];
        assert( jblk0<(int)iblk );
        const T* vij = &vcrs[ijcrs*9];
        const T valj0 = vec[jblk0*3+0];
        const T valj1 = vec[jblk0*3+1];
        const T valj2 = vec[jblk0*3+2];
        pTmpVec[0] -= vij[0]*valj0+vij[1]*valj1+vij[2]*valj2;
        pTmpVec[1] -= vij[3]*valj0+vij[4]*valj1+vij[5]*valj2;
        pTmpVec[2] -= vij[6]*valj0+vij[7]*valj1+vij[8]*valj2;
      }
      const T* vii = &vdia[iblk*9];
      vec[iblk*3+0] = vii[0]*pTmpVec[0]+vii[1]*pTmpVec[1]+vii[2]*pTmpVec[2];
      vec[iblk*3+1] = vii[3]*pTmpVec[0]+vii[4]*pTmpVec[1]+vii[5]*pTmpVec[2];
      vec[iblk*3+2] = vii[6]*pTmpVec[0]+vii[7]*pTmpVec[1]+vii[8]*pTmpVec[2];
    }
  }
  else if (len==4){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    const T* vdia = mat.valDia.data();
    // ------------
    T pTmpVec[4];
    for (unsigned int iblk = 0; iblk<nblk; iblk++){
      pTmpVec[0] = vec[iblk*4+0];
      pTmpVec[1] = vec[iblk*4+1];
      pTmpVec[2] = vec[iblk*4+2];
      pTmpVec[3] = vec[iblk*4+3];
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = m_diaInd[iblk];
      for (unsigned int ijcrs = icrs0; ijcrs<icrs1; ijcrs++){
        assert(ijcrs<mat.rowPtr.size());
        const int jblk0 = rowptr[ijcrs];
        assert(jblk0<(int)iblk);
        const T* vij = &vcrs[ijcrs*16];
        const T valj0 = vec[jblk0*4+0];
        const T valj1 = vec[jblk0*4+1];
        const T valj2 = vec[jblk0*4+2];
        const T valj3 = vec[jblk0*4+3];
        pTmpVec[0] -= vij[ 0]*valj0+vij[ 1]*valj1+vij[ 2]*valj2+vij[ 3]*valj3;
        pTmpVec[1] -= vij[ 4]*valj0+vij[ 5]*valj1+vij[ 6]*valj2+vij[ 7]*valj3;
        pTmpVec[2] -= vij[ 8]*valj0+vij[ 9]*valj1+vij[10]*valj2+vij[11]*valj3;
        pTmpVec[3] -= vij[12]*valj0+vij[13]*valj1+vij[14]*valj2+vij[15]*valj3;
      }
      const T* vii = &vdia[iblk*16];
      vec[iblk*4+0] = vii[ 0]*pTmpVec[0]+vii[ 1]*pTmpVec[1]+vii[ 2]*pTmpVec[2]+vii[ 3]*pTmpVec[3];
      vec[iblk*4+1] = vii[ 4]*pTmpVec[0]+vii[ 5]*pTmpVec[1]+vii[ 6]*pTmpVec[2]+vii[ 7]*pTmpVec[3];
      vec[iblk*4+2] = vii[ 8]*pTmpVec[0]+vii[ 9]*pTmpVec[1]+vii[10]*pTmpVec[2]+vii[11]*pTmpVec[3];
      vec[iblk*4+3] = vii[12]*pTmpVec[0]+vii[13]*pTmpVec[1]+vii[14]*pTmpVec[2]+vii[15]*pTmpVec[3];
    }
  }
  else{
    const int blksize = len*len;
    std::vector<T> pTmpVec(len);
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      for(unsigned int idof=0;idof<len;idof++){
        pTmpVec[idof] = vec[iblk*len+idof];
      }
      for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<m_diaInd[iblk];ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = mat.rowPtr[ijcrs];
        assert( jblk0<(int)iblk );
        const T* vij = &mat.valCrs[ijcrs*blksize];
        for(unsigned int idof=0;idof<len;idof++){
          for(unsigned int jdof=0;jdof<len;jdof++){
            pTmpVec[idof] -= vij[idof*len+jdof]*vec[jblk0*len+jdof];
          }
        }
      }
      const T* vii = &mat.valDia[iblk*blksize];
      for(unsigned int idof=0;idof<len;idof++){
        T dtmp1 = 0.0;
        for(unsigned int jdof=0;jdof<len;jdof++){
          dtmp1 += vii[idof*len+jdof]*pTmpVec[jdof];
        }
        vec[iblk*len+idof] = dtmp1;
      }
    }
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CPreconditionerILU<double>::ForwardSubstitution(
    double* vec ) const;
template void delfem2::CPreconditionerILU<std::complex<double>>::ForwardSubstitution(
    std::complex<double>* vec ) const;
#endif

// ------------------------------------------------------------------

template <typename T>
void delfem2::CPreconditionerILU<T>::ForwardSubstitutionDegenerate(
    T* vec,
    const unsigned int N ) const
{
  assert( mat.nrowdim == 1 );
  const unsigned int nblk = mat.nrowblk;
  
  const unsigned int* colind = mat.colInd.data();
  const unsigned int* rowptr = mat.rowPtr.data();
  const T* vcrs = mat.valCrs.data();
  const T* vdia = mat.valDia.data();
  // -------------------------
  std::vector<T> buff(N);
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    memcpy(buff.data(), vec+iblk*N, N*sizeof(T));
    for(unsigned int ijcrs=colind[iblk];ijcrs<m_diaInd[iblk];ijcrs++){
      assert( ijcrs<mat.rowPtr.size() );
      const int jblk0 = rowptr[ijcrs];
      assert( jblk0<(int)iblk );
      for(unsigned int i=0;i<N;++i){
        buff[i] -= vcrs[ijcrs]*vec[jblk0*N+i];
      }
    }
    for(unsigned int i=0;i<N;++i){
      vec[iblk*N+i] = vdia[iblk]*buff[i];
    }
  }
}


// --------------------------------------------------

template <typename T>
void delfem2::CPreconditionerILU<T>::BackwardSubstitution
( T* vec ) const
{
  const unsigned int len = mat.nrowdim;
  const int nblk = mat.nrowblk;
  
  if( len == 1 ){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    // -------------------------------
    for(int iblk=nblk-1;iblk>=0;iblk--){
      assert( (int)iblk < nblk );
      T lvec_i = vec[iblk];
      for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = rowptr[ijcrs];
        assert( jblk0>(int)iblk && jblk0<nblk );
        lvec_i -= vcrs[ijcrs]*vec[jblk0];
      }
      vec[iblk] = lvec_i;
    }
  }
  else if( len == 2 ){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    // ----------------------------
    T pTmpVec[2];
    for(int iblk=nblk-1;iblk>=0;iblk--){
      assert( (int)iblk < nblk );
      pTmpVec[0] = vec[iblk*2+0];
      pTmpVec[1] = vec[iblk*2+1];
      const unsigned int icrs0 = m_diaInd[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = rowptr[ijcrs];
        assert( jblk0>(int)iblk && jblk0<nblk );
        const T* vij = &vcrs[ijcrs*4];
        const T valj0 = vec[jblk0*2+0];
        const T valj1 = vec[jblk0*2+1];
        pTmpVec[0] -= vij[0]*valj0+vij[1]*valj1;
        pTmpVec[1] -= vij[2]*valj0+vij[3]*valj1;
      }
      vec[iblk*2+0] = pTmpVec[0];
      vec[iblk*2+1] = pTmpVec[1];
    }
  }
  else if( len == 3 ){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    // --------------------
    T pTmpVec[3];
    for(int iblk=nblk-1;iblk>=0;iblk--){
      assert( (int)iblk < nblk );
      pTmpVec[0] = vec[iblk*3+0];
      pTmpVec[1] = vec[iblk*3+1];
      pTmpVec[2] = vec[iblk*3+2];
      const int icrs0 = m_diaInd[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = rowptr[ijcrs];
        assert( jblk0>(int)iblk && jblk0<nblk );
        const T* vij = &vcrs[ijcrs*9];
        const T valj0 = vec[jblk0*3+0];
        const T valj1 = vec[jblk0*3+1];
        const T valj2 = vec[jblk0*3+2];
        pTmpVec[0] -= vij[0]*valj0+vij[1]*valj1+vij[2]*valj2;
        pTmpVec[1] -= vij[3]*valj0+vij[4]*valj1+vij[5]*valj2;
        pTmpVec[2] -= vij[6]*valj0+vij[7]*valj1+vij[8]*valj2;
      }
      vec[iblk*3+0] = pTmpVec[0];
      vec[iblk*3+1] = pTmpVec[1];
      vec[iblk*3+2] = pTmpVec[2];
    }
  }
  else if (len==4){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    // -----------------------------
    T pTmpVec[4];
    for (int iblk = nblk-1; iblk>=0; iblk--){
      assert((int)iblk < nblk);
      pTmpVec[0] = vec[iblk*4+0];
      pTmpVec[1] = vec[iblk*4+1];
      pTmpVec[2] = vec[iblk*4+2];
      pTmpVec[3] = vec[iblk*4+3];
      const int icrs0 = m_diaInd[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for (unsigned int ijcrs = icrs0; ijcrs<icrs1; ijcrs++){
        assert(ijcrs<mat.rowPtr.size());
        const int jblk0 = rowptr[ijcrs];
        assert(jblk0>(int)iblk && jblk0<nblk);
        const T* vij = &vcrs[ijcrs*16];
        const T valj0 = vec[jblk0*4+0];
        const T valj1 = vec[jblk0*4+1];
        const T valj2 = vec[jblk0*4+2];
        const T valj3 = vec[jblk0*4+3];
        pTmpVec[0] -= vij[ 0]*valj0+vij[ 1]*valj1+vij[ 2]*valj2+vij[ 3]*valj3;
        pTmpVec[1] -= vij[ 4]*valj0+vij[ 5]*valj1+vij[ 6]*valj2+vij[ 7]*valj3;
        pTmpVec[2] -= vij[ 8]*valj0+vij[ 9]*valj1+vij[10]*valj2+vij[11]*valj3;
        pTmpVec[3] -= vij[12]*valj0+vij[13]*valj1+vij[14]*valj2+vij[15]*valj3;
      }
      vec[iblk*4+0] = pTmpVec[0];
      vec[iblk*4+1] = pTmpVec[1];
      vec[iblk*4+2] = pTmpVec[2];
      vec[iblk*4+3] = pTmpVec[3];
    }
  }
  else{
    const int blksize = len*len;
    std::vector<T> pTmpVec(len);
    for(int iblk=nblk-1;iblk>=0;iblk--){
      assert( (int)iblk < nblk );
      for(unsigned int idof=0;idof<len;idof++){
        pTmpVec[idof] = vec[iblk*len+idof];
      }
      for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = mat.rowPtr[ijcrs];
        assert( jblk0>(int)iblk && jblk0<nblk );
        const T* vij = &mat.valCrs[ijcrs*blksize];
        for(unsigned int idof=0;idof<len;idof++){
          for(unsigned int jdof=0;jdof<len;jdof++){
            pTmpVec[idof] -= vij[idof*len+jdof]*vec[jblk0*len+jdof];
          }
        }
      }
      for(unsigned int idof=0;idof<len;idof++){
        vec[iblk*len+idof] = pTmpVec[idof];
      }
    }
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CPreconditionerILU<double>::BackwardSubstitution(
    double* vec ) const;
template void delfem2::CPreconditionerILU<std::complex<double>>::BackwardSubstitution(
    std::complex<double>* vec ) const;
#endif

// -----------------------------------------------------------------


template <typename T>
void delfem2::CPreconditionerILU<T>::BackwardSubstitutionDegenerate
 ( T* vec, unsigned int N ) const
{
  assert( mat.nrowdim == 1 );
  const int nblk = mat.nrowblk;
  
  const unsigned int* colind = mat.colInd.data();
  const unsigned int* rowptr = mat.rowPtr.data();
  const T* vcrs = mat.valCrs.data();
  
  std::vector<T> buff(N);
  for(int iblk=nblk-1;iblk>=0;iblk--){
    assert( (int)iblk < nblk );
    memcpy(buff.data(), vec+iblk*N, N*sizeof(T));
    for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<colind[iblk+1];ijcrs++){
      assert( ijcrs<mat.rowPtr.size() );
      const int jblk0 = rowptr[ijcrs];
      assert( jblk0>(int)iblk && jblk0<nblk );
      for(unsigned int i=0;i<N;++i){
        buff[i] -= vcrs[ijcrs]*vec[jblk0*N+i];
      }
    }
    for(unsigned int i=0;i<N;++i){
      vec[iblk*N+i] = buff[i];
    }
  }
}
  
// --------------------------------------------------------------

// if(lev_fill == -1){ take all the fills }
template <typename T>
void delfem2::CPreconditionerILU<T>::Initialize_ILUk
(const CMatrixSparse<T>& m,
 int lev_fill)
{
  
  if (lev_fill==0){
    this->Initialize_ILU0(m);
    return;
  }
    
  std::vector<ilu::CRowLev> aRowLev;
  aRowLev.reserve(m.rowPtr.size()*4);
  
  assert(m.nrowblk==m.ncolblk);
  const unsigned int nblk = m.nrowblk;
  assert(m.nrowdim==m.ncoldim);
  const unsigned int len = m.nrowdim;
  
  assert(!m.valDia.empty());
  mat.Initialize(nblk, len, true);
  
  m_diaInd.resize(nblk);
  
  for(unsigned int iblk=0; iblk<nblk; ++iblk){
    std::vector<ilu::CRowLevNext> listNonzero;
    {  // copy row pattern of input matrix into listNonzero
      listNonzero.resize(m.colInd[iblk+1]-m.colInd[iblk]);
      int inz = 0;
      for (unsigned int ijcrs = m.colInd[iblk]; ijcrs<m.colInd[iblk+1]; ijcrs++){
        assert(ijcrs<m.rowPtr.size());
        const unsigned int jblk0 = m.rowPtr[ijcrs];
        assert(jblk0<nblk);
        listNonzero[inz].row = jblk0;
        listNonzero[inz].lev = 0;
        listNonzero[inz].next = inz+1;
        inz++;
      }
      listNonzero[inz-1].next = -1;
    }
    
    int knz_cur = 0;
    for (;;){
      const unsigned int kblk0 = listNonzero[knz_cur].row;
      assert(kblk0<nblk);
      const unsigned int ik_lev0 = listNonzero[knz_cur].lev;
      if ((int)ik_lev0+1>lev_fill && lev_fill!=-1){
        knz_cur = listNonzero[knz_cur].next;
        if (knz_cur==-1) break;
        continue;
      }
      if (kblk0>=iblk) break;
      
      int jnz_cur = knz_cur;
      for (unsigned int kjcrs = m_diaInd[kblk0]; kjcrs<mat.colInd[kblk0+1]; kjcrs++){
        const unsigned int kj_lev0 = aRowLev[kjcrs].lev;
        if ((int)kj_lev0+1>lev_fill && lev_fill!=-1) continue;
        const unsigned int jblk0 = aRowLev[kjcrs].row;
        assert(jblk0>kblk0 && jblk0<nblk);
        assert(listNonzero[jnz_cur].row < jblk0);
        if (jblk0==iblk) continue; // already filled-in on the diagonal
        
        // check if this is fill in
        bool is_fill_in = false;
        for (;;){
          const int jnz_nex = listNonzero[jnz_cur].next;
          assert((jnz_nex>=0&&jnz_nex<(int)nblk)||jnz_nex==-1);
          if (jnz_nex==-1){ is_fill_in = true; break; }
          if (listNonzero[jnz_nex].row>jblk0){ is_fill_in = true; break; }
          if (listNonzero[jnz_nex].row==jblk0){ break; }
          assert(listNonzero[jnz_nex].row < jblk0);
          jnz_cur = jnz_nex;
        }
        if (!is_fill_in){ continue; }
        
        // pick up fill in
        const unsigned int max_lev0 = (ik_lev0 > kj_lev0) ? ik_lev0 : kj_lev0;
        const unsigned  int inz_last = listNonzero.size();
        listNonzero.resize(listNonzero.size()+1);
        listNonzero[inz_last].row = jblk0;
        listNonzero[inz_last].lev = max_lev0+1;
        listNonzero[inz_last].next = listNonzero[jnz_cur].next;
        listNonzero[jnz_cur].next = inz_last;
        jnz_cur = inz_last;
      }
      knz_cur = listNonzero[knz_cur].next;
      assert((knz_cur>=0&&knz_cur<(int)nblk)||knz_cur==-1);
      if (knz_cur==-1) break;
    }
    
    // -------------------------------------
    
    {
      aRowLev.resize(mat.colInd[iblk]+listNonzero.size());
      unsigned int icrs0 = mat.colInd[iblk];
      for (int inz = 0; inz!=-1; inz = listNonzero[inz].next){
        const unsigned int jblk = listNonzero[inz].row;
        const unsigned int jlev = listNonzero[inz].lev;
        assert(jblk<nblk);
        assert(jblk!=iblk);
        aRowLev[icrs0].row = jblk;
        aRowLev[icrs0].lev = jlev;
        icrs0++;
      }
      
      mat.colInd[iblk+1] = icrs0;
      mat.rowPtr.resize(icrs0);
      m_diaInd[iblk] = icrs0;
      for (unsigned int ijcrs = mat.colInd[iblk]; ijcrs<mat.colInd[iblk+1]; ijcrs++){
        const unsigned int jblk0 = aRowLev[ijcrs].row;
        if (jblk0 > iblk){
          m_diaInd[iblk] = ijcrs;
          break;
        }
      }
    }
  }
  
  {
    const unsigned int ncrs = mat.rowPtr.size();
    //    std::cout << aRowLev.size() << " " << ncrs << std::endl;
    assert(aRowLev.size()==ncrs);
    mat.rowPtr.resize(ncrs);
    for (unsigned int icrs = 0; icrs<ncrs; ++icrs){
      mat.rowPtr[icrs] = aRowLev[icrs].row;
    }
    const int blksize = len*len;
    mat.valCrs.resize(ncrs*blksize);
    assert(!mat.valDia.empty());
    mat.valDia = m.valDia;
    //    std::cout<<"ncrs: "<<ncrs<<" "<<m.rowPtr.size()<<std::endl;
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CPreconditionerILU<double>::Initialize_ILUk(
    const CMatrixSparse<double>& m, int lev_fill);
template void delfem2::CPreconditionerILU<std::complex<double>>::Initialize_ILUk(
    const CMatrixSparse<std::complex<double>>& m, int lev_fill);
#endif


template <typename T>
void delfem2::CPreconditionerILU<T>::SetValueILU
 (const CMatrixSparse<T>& rhs)
{
  const unsigned int nblk = mat.nrowblk;
  const unsigned int len = mat.nrowdim;
  assert( rhs.nrowblk == nblk );
  assert( rhs.ncolblk == nblk );
  assert( rhs.nrowdim == len );
  assert( rhs.ncoldim == len );
  const unsigned int blksize = len*len;
  std::vector<int> row2crs(nblk,-1);
  {
    const unsigned int n = mat.rowPtr.size()*len*len;
    for(unsigned int i=0;i<n;++i){ mat.valCrs[i] = 0.0; }
  }
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<mat.rowPtr.size() );
      const unsigned int jblk0 = mat.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = ijcrs;
    }
    for(unsigned int ijcrs=rhs.colInd[iblk];ijcrs<rhs.colInd[iblk+1];ijcrs++){
      assert( ijcrs<rhs.rowPtr.size() );
      const unsigned int jblk0 = rhs.rowPtr[ijcrs];
      assert( jblk0<nblk );
      const int ijcrs0 = row2crs[jblk0];
      if( ijcrs0 == -1 ) continue;
      const T* pval_in = &rhs.valCrs[ijcrs*blksize];
      T* pval_out = &mat.valCrs[ijcrs0*blksize];
      for(unsigned int i=0;i<blksize;i++){ *(pval_out+i) = *(pval_in+i); }
    }
    for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<mat.rowPtr.size() );
      const unsigned int jblk0 = mat.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = -1;
    }
  }
  for(unsigned int i=0;i<nblk*blksize;i++){ mat.valDia[i] = rhs.valDia[i]; }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CPreconditionerILU<double>::SetValueILU(
    const CMatrixSparse<double>& rhs);
template void delfem2::CPreconditionerILU<std::complex<double>>::SetValueILU(
    const CMatrixSparse<std::complex<double>>& rhs);
#endif


template <typename T>
void delfem2::CPreconditionerILU<T>::Initialize_ILU0
 (const CMatrixSparse<T>& m)
{
  this->mat = m;
  const int nblk = m.nrowblk;
  m_diaInd.resize(nblk);
  for(int iblk=0;iblk<nblk;iblk++){
    m_diaInd[iblk] = mat.colInd[iblk+1];
    for(unsigned int icrs=mat.colInd[iblk];icrs<mat.colInd[iblk+1];icrs++){
      assert( icrs < mat.rowPtr.size() );
      const int jblk0 = mat.rowPtr[icrs];
      assert( jblk0 < nblk );
      if( jblk0 > iblk ){
        m_diaInd[iblk] = icrs;
        break;
      }
    }
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CPreconditionerILU<double>::Initialize_ILU0(
    const CMatrixSparse<double>& m);
template void delfem2::CPreconditionerILU<std::complex<double>>::Initialize_ILU0(
    const CMatrixSparse<std::complex<double>>& m);
#endif

