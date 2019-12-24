/**
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <complex>

#include "delfem2/ilu_mats.h"

typedef std::complex<double> COMPLEX;
namespace dfm2 = delfem2;

// ----------------------------------------------------

static void CalcMatPr(double* out, const double* d, double* tmp,
                      const unsigned int ni, const unsigned int nj )
{
	unsigned int i,j,k;
	for(i=0;i<ni;i++){
		for(j=0;j<nj;j++){
			tmp[i*nj+j] = 0.0;
			for(k=0;k<ni;k++){
				tmp[i*nj+j] += d[i*ni+k]*out[k*nj+j];
			}
		}
	}
	for(i=0;i<ni*nj;i++){
		out[i] = tmp[i];
	}
}

static void CalcSubMatPr(double* out, const double* a, const double* b,
                         const int ni, const int nk, const int nj )
{
	int i,j,k;
	for(i=0;i<ni;i++){
		for(j=0;j<nj;j++){
			for(k=0;k<nk;k++){
				out[i*nj+j] -= a[i*nk+k]*b[k*nj+j];
			}
		}
	}
}


static void CalcInvMat(
    double* a,
    const unsigned int n,
    int& info )
{
	double tmp1;
  
	info = 0;
	unsigned int i,j,k;
	for(i=0;i<n;i++){
		if( fabs(a[i*n+i]) < 1.0e-30 ){
			info = 1;
			return;
		}
		if( a[i*n+i] < 0.0 ){
			info--;
		}
		tmp1 = 1.0 / a[i*n+i];
		a[i*n+i] = 1.0;
		for(k=0;k<n;k++){
			a[i*n+k] *= tmp1;
		}
		for(j=0;j<n;j++){
			if( j!=i ){
				tmp1 = a[j*n+i];
				a[j*n+i] = 0.0;
				for(k=0;k<n;k++){
					a[j*n+k] -= tmp1*a[i*n+k];
				}
			}
		}
	}
}

// t is a tmporary buffer size of 9
static inline void CalcInvMat3(double a[], double t[] )
{
	const double det =
  + a[0]*a[4]*a[8] + a[3]*a[7]*a[2] + a[6]*a[1]*a[5]
  - a[0]*a[7]*a[5] - a[6]*a[4]*a[2] - a[3]*a[1]*a[8];
	const double inv_det = 1.0/det;
  
  for(int i=0;i<9;i++){ t[i] = a[i]; }
  
	a[0] = inv_det*(t[4]*t[8]-t[5]*t[7]);
	a[1] = inv_det*(t[2]*t[7]-t[1]*t[8]);
	a[2] = inv_det*(t[1]*t[5]-t[2]*t[4]);
  
	a[3] = inv_det*(t[5]*t[6]-t[3]*t[8]);
	a[4] = inv_det*(t[0]*t[8]-t[2]*t[6]);
	a[5] = inv_det*(t[2]*t[3]-t[0]*t[5]);
  
	a[6] = inv_det*(t[3]*t[7]-t[4]*t[6]);
	a[7] = inv_det*(t[1]*t[6]-t[0]*t[7]);
	a[8] = inv_det*(t[0]*t[4]-t[1]*t[3]);
}

/* --------------------------------------------------------------------- */


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

// numerical factorization
template <>
bool delfem2::CPreconditionerILU<double>::DoILUDecomp()
{
  const int nmax_sing = 10;
	int icnt_sing = 0;
  
	const unsigned int len = mat.len_col;
  const unsigned int nblk = mat.nblk_col;
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
					CalcInvMat3(vii,tmpBlk);
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
          CalcSubMatPr(vij,vik,vkj, len,len,len);
				}
			}
			{
				double* vii = &vdia[iblk*blksize];
				int info = 0;
				CalcInvMat(vii,len,info);
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
        CalcMatPr(vij,pVal_ii,pTmpBlk,  len,len);
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
bool delfem2::CPreconditionerILU<COMPLEX>::DoILUDecomp()
{
//  const int nmax_sing = 10;
//  int icnt_sing = 0;
  
  const unsigned int len = mat.len_col;
  const unsigned int nblk = mat.nblk_col;
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


// -----------------------------------------------

template <>
std::vector<double> delfem2::Solve_PCG
(double* r_vec,
 double* x_vec,
 double conv_ratio_tol,
 unsigned int max_nitr,
 const CMatrixSparse<double>& mat,
 const delfem2::CPreconditionerILU<double>& ilu)
{
  const unsigned int ndof = mat.nblk_col*mat.len_col;
  std::vector<double> aResHistry;
  
  for(unsigned int i=0;i<ndof;i++){ x_vec[i] = 0; }    // {x} = 0
  
	double inv_sqnorm_res0;
	{
		const double sqnorm_res0 = DotX(r_vec,r_vec,ndof);
    aResHistry.push_back(sqrt(sqnorm_res0));
		if( sqnorm_res0 < 1.0e-30 ){ return aResHistry; }
		inv_sqnorm_res0 = 1.0 / sqnorm_res0;
	}
  
  // {Pr} = [P]{r}
  std::vector<double> Pr_vec(r_vec,r_vec+ndof);
  ilu.Solve(Pr_vec);
  // {p} = {Pr}
  std::vector<double> p_vec = Pr_vec;
  // rPr = ({r},{Pr})
	double rPr = DotX(r_vec,Pr_vec.data(),ndof);
	for(unsigned int iitr=0;iitr<max_nitr;iitr++){
		{
      std::vector<double>& Ap_vec = Pr_vec;      
      // {Ap} = [A]{p}
			mat.MatVec(1.0,p_vec.data(),0.0,Ap_vec.data());
      // alpha = ({r},{Pr})/({p},{Ap})
			const double pAp = Dot(p_vec,Ap_vec);
			double alpha = rPr / pAp;
      AXPY(-alpha,Ap_vec.data(),r_vec, ndof);       // {r} = -alpha*{Ap} + {r}
      AXPY(+alpha,p_vec.data(), x_vec, ndof);       // {x} = +alpha*{p } + {x}
    }
		{	// Converge Judgement
			double sqnorm_res = DotX(r_vec,r_vec,ndof);
      aResHistry.push_back(sqrt(sqnorm_res));
      double conv_ratio = sqrt(sqnorm_res*inv_sqnorm_res0);
      if( conv_ratio < conv_ratio_tol ){ return aResHistry; }
		}
		{	// calc beta
      // {Pr} = [P]{r}
      for(unsigned int i=0;i<ndof;i++){ Pr_vec[i] = r_vec[i]; }
			ilu.Solve(Pr_vec);
      // rPr1 = ({r},{Pr})
			const double rPr1 = DotX(r_vec,Pr_vec.data(),ndof);
      // beta = rPr1/rPr
			double beta = rPr1/rPr;
			rPr = rPr1;
      // {p} = {Pr} + beta*{p}
      for(unsigned int i=0;i<ndof;i++){ p_vec[i] = Pr_vec[i] + beta*p_vec[i]; }
    }
	}
  {
    // Converge Judgement
    double sq_norm_res = DotX(r_vec,r_vec,ndof);
    aResHistry.push_back(sqrt(sq_norm_res));
  }
  return aResHistry;
}


template <>
std::vector<double> delfem2::Solve_PCG
(COMPLEX* r_vec,
 COMPLEX* x_vec,
 double conv_ratio_tol,
 unsigned int max_nitr,
 const CMatrixSparse<COMPLEX>& mat,
 const delfem2::CPreconditionerILU<COMPLEX>& ilu)
{
  const unsigned int ndof = mat.nblk_col*mat.len_col;
  std::vector<double> aResHistry;
  
  for(unsigned int i=0;i<ndof;i++){ x_vec[i] = COMPLEX(0.0,0.0); }    // {x} = 0
  
  double inv_sqnorm_res0;
  {
    const double sqnorm_res0 = DotX(r_vec,r_vec,ndof).real();
    aResHistry.push_back(sqnorm_res0);
    if( sqnorm_res0 < 1.0e-30 ){ return aResHistry; }
    inv_sqnorm_res0 = 1.0 / sqnorm_res0;
  }
  
  // {Pr} = [P]{r}
  std::vector<COMPLEX> Pr_vec(r_vec,r_vec+ndof);
  ilu.Solve(Pr_vec);
  // {p} = {Pr}
  std::vector<COMPLEX> p_vec = Pr_vec;
  // rPr = ({r},{Pr})
  COMPLEX rPr = DotX(r_vec,Pr_vec.data(),ndof);
  for(unsigned int iitr=0;iitr<max_nitr;iitr++){
    {
      std::vector<COMPLEX>& Ap_vec = Pr_vec;
      // {Ap} = [A]{p}
      mat.MatVec(1.0,p_vec.data(),0.0,Ap_vec.data());
      // alpha = ({r},{Pr})/({p},{Ap})
      const double pAp = Dot(p_vec,Ap_vec).real();
      COMPLEX alpha = rPr / pAp;
      AXPY(-alpha,Ap_vec.data(),r_vec, ndof);       // {r} = -alpha*{Ap} + {r}
      AXPY(+alpha,p_vec.data(), x_vec, ndof);       // {x} = +alpha*{p } + {x}
    }
    {  // Converge Judgement
      double sqnorm_res = DotX(r_vec,r_vec,ndof).real();
      double conv_ratio = sqrt(sqnorm_res*inv_sqnorm_res0);
      aResHistry.push_back(conv_ratio);
      if( conv_ratio < conv_ratio_tol ){ return aResHistry; }
    }
    {  // calc beta
      // {Pr} = [P]{r}
      for(unsigned int i=0;i<ndof;i++){ Pr_vec[i] = r_vec[i]; }
      ilu.Solve(Pr_vec);
      // rPr1 = ({r},{Pr})
      const COMPLEX rPr1 = DotX(r_vec,Pr_vec.data(),ndof);
      // beta = rPr1/rPr
      COMPLEX beta = rPr1/rPr;
      rPr = rPr1;
      // {p} = {Pr} + beta*{p}
      for(unsigned int i=0;i<ndof;i++){ p_vec[i] = Pr_vec[i] + beta*p_vec[i]; }
    }
  }
  {
    // Converge Judgement
    double sq_norm_res = DotX(r_vec,r_vec,ndof).real();
    aResHistry.push_back(sqrt(sq_norm_res));
  }
  return aResHistry;
}

template <>
std::vector<double> delfem2::Solve_PBiCGStab
(double* r_vec,
 double* x_vec,
 double conv_ratio_tol,
 unsigned int max_niter,
 const CMatrixSparse<double>& mat,
 const delfem2::CPreconditionerILU<double>& ilu)
{
  assert( !mat.valDia.empty() );
  assert( mat.nblk_col == mat.nblk_row );
  assert( mat.len_col == mat.len_row );
  const unsigned int ndof = mat.nblk_col*mat.len_col;
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
    ilu.Solve(Mp_vec);
    // calc (r,r0*)
    const double r_r2 = DotX(r_vec,r0_vec.data(),ndof);
    // calc {AMp_vec} = [A]*{Mp_vec}
    mat.MatVec(1.0, Mp_vec.data(), 0.0, AMp_vec.data());
    // calc alpha
    const double alpha = r_r2 / Dot(AMp_vec,r0_vec);
    // calc s_vector
    s_vec.assign(r_vec,r_vec+ndof);
    AXPY(-alpha,AMp_vec,s_vec);
    // {Ms_vec} = [M^-1]*{s}
    Ms_vec = s_vec;
    ilu.Solve(Ms_vec);
    // calc {AMs_vec} = [A]*{Ms_vec}
    mat.MatVec(1.0,Ms_vec.data(),0.0,AMs_vec.data());
    double omega;
    {	// calc omega
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
    {	// calc beta
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



template <>
std::vector<double> delfem2::Solve_PBiCGStab
(COMPLEX* r_vec,
 COMPLEX* x_vec,
 double conv_ratio_tol,
 unsigned int max_niter,
 const CMatrixSparse<COMPLEX>& mat,
 const delfem2::CPreconditionerILU<COMPLEX>& ilu)
{
  assert( !mat.valDia.empty() );
  assert( mat.nblk_col == mat.nblk_row );
  assert( mat.len_col == mat.len_row );
  const unsigned int ndof = mat.nblk_col*mat.len_col;
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
  
  std::vector<COMPLEX> s_vec(ndof);
  std::vector<COMPLEX> Ms_vec(ndof);
  std::vector<COMPLEX> AMs_vec(ndof);
  std::vector<COMPLEX> Mp_vec(ndof);
  std::vector<COMPLEX> AMp_vec(ndof);
  
  const std::vector<COMPLEX> r0_vec(r_vec,r_vec+ndof);   // {r2} = {r}
  std::vector<COMPLEX> p_vec(r_vec,r_vec+ndof);  // {p} = {r}
  
  // calc (r,r0*)
  COMPLEX r_r0 = DotX(r_vec,r0_vec.data(),ndof);
  
  for(unsigned int itr=0;itr<max_niter;itr++){
    // {Mp_vec} = [M^-1]*{p}
    Mp_vec.assign(p_vec.begin(),p_vec.end());
    ilu.Solve(Mp_vec);
    // calc {AMp_vec} = [A]*{Mp_vec}
    mat.MatVec(COMPLEX(1,0), Mp_vec.data(), COMPLEX(0,0), AMp_vec.data());
    // calc alpha
    const COMPLEX alpha = r_r0 / Dot(AMp_vec,r0_vec);
    // calc s_vector
    s_vec.assign(r_vec,r_vec+ndof);
    AXPY(-alpha,AMp_vec,s_vec);
    // {Ms_vec} = [M^-1]*{s}
    Ms_vec.assign(s_vec.begin(),s_vec.end());
    ilu.Solve(Ms_vec);
    // calc {AMs_vec} = [A]*{Ms_vec}
    mat.MatVec(COMPLEX(1,0),Ms_vec.data(), COMPLEX(0,0), AMs_vec.data());
    const COMPLEX omega = Dot(s_vec,AMs_vec) / Dot(AMs_vec,AMs_vec).real();
    for(unsigned int i=0;i<ndof;++i){ x_vec[i] = x_vec[i]+alpha*Mp_vec[i]+omega*Ms_vec[i]; }
    for(unsigned int i=0;i<ndof;++i){ r_vec[i] = s_vec[i]-omega*AMs_vec[i]; }
    {
      const double sq_norm_res = DotX(r_vec,r_vec,ndof).real();
      const double conv_ratio = sqrt(sq_norm_res * sq_inv_norm_res_ini);
      aResHistry.push_back( conv_ratio );
      if( conv_ratio < conv_ratio_tol ){ return aResHistry; }
    }
    COMPLEX beta;
    {  // calc beta
      const COMPLEX tmp1 = DotX(r_vec,r0_vec.data(),ndof);
      beta = (tmp1*alpha)/(r_r0*omega);
      r_r0 = tmp1;
    }
    // update p_vector
    for(unsigned int i=0;i<ndof;++i){ p_vec[i] = r_vec[i]+beta*(p_vec[i]-omega*AMp_vec[i]); }
  }
  
  return aResHistry;
}


std::vector<double> dfm2::Solve_PCOCG
(COMPLEX* r_vec,
 COMPLEX* x_vec,
 double conv_ratio_tol,
 unsigned int max_niter,
 const CMatrixSparse<COMPLEX>& mat,
 const delfem2::CPreconditionerILU<COMPLEX>& ilu)
{
  assert( !mat.valDia.empty() );
  assert( mat.nblk_col == mat.nblk_row );
  assert( mat.len_col == mat.len_row );
  const unsigned int ndof = mat.nblk_col*mat.len_col;
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
  ilu.Solve(w_vec);
  
  std::vector<COMPLEX> p_vec = w_vec;  // {p} = {w}
  COMPLEX r_w = MultSumX(r_vec,w_vec.data(),ndof);
  
  for(unsigned int itr=0;itr<max_niter;itr++){
    mat.MatVec(COMPLEX(1,0), p_vec.data(), COMPLEX(0,0), Ap_vec.data());
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
    ilu.Solve(w_vec);
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


