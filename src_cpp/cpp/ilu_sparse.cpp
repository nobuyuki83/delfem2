/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <assert.h>
#include <math.h>
#include <vector>

#include "delfem2/ilu_sparse.h"


static void CalcMatPr(double* out, const double* d, double* tmp,
                      const int ni, const int nj )
{
	int i,j,k;
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


static void CalcInvMat(double* a, const int n, int& info )
{
	double tmp1;
  
	info = 0;
	int i,j,k;
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


CPreconditionerILU::CPreconditionerILU()
{
  std::cout << "CPreconditionerILU -- construct" << std::endl;
}


CPreconditionerILU::CPreconditionerILU(const CPreconditionerILU& p)
{
  std::cout << "CPreconditionerILU -- construct copy" << std::endl;
  this->mat = p.mat; // deep copy
  const int nblk = this->mat.nblk_col;
  this->m_diaInd.resize(nblk);
  for(int iblk=0;iblk<nblk;++iblk){
    this->m_diaInd[iblk] = p.m_diaInd[iblk];
  }
}



CPreconditionerILU::~CPreconditionerILU()
{
  std::cout << "CPreconditionerILU -- destroy" << std::endl;
  m_diaInd.clear();
  std::cout << "CPreconditionerILU -- destroy end" << std::endl;
}


void CPreconditionerILU::Initialize_ILU0
(const CMatrixSparse& m)
{
  this->mat = m;
  const int nblk = m.nblk_col;

  m_diaInd.resize(nblk);
  
  for(int iblk=0;iblk<nblk;iblk++){
    m_diaInd[iblk] = mat.colInd[iblk+1];
    for(unsigned int icrs=mat.colInd[iblk];icrs<mat.colInd[iblk+1];icrs++){
      assert( icrs < mat.ncrs );
      const int jblk0 = mat.rowPtr[icrs];
      assert( jblk0 < nblk );
      if( jblk0 > iblk ){
        m_diaInd[iblk] = icrs;
        break;
      }
    }
  }
}


class CRowLev{
public:
  CRowLev() :row(0), lev(0){}
  CRowLev(int row, int lev) : row(row), lev(lev){}
  bool operator < (const CRowLev& rhs) const{
    if (row!=rhs.row) return row < rhs.row;
    return lev < rhs.lev;
  }
public:
  int row;
  int lev;
};

class CRowLevNext{
public:
  int row;
  int lev;
  int next;
};

// if(lev_fill == -1){ take all the fills }
void CPreconditionerILU::Initialize_ILUk
(const CMatrixSparse& m,
int lev_fill)
{

  if (lev_fill==0){
   this->Initialize_ILU0(m);
   return;
  }

  std::vector<CRowLev> aRowLev;
  aRowLev.reserve(m.ncrs*4);


  assert(m.nblk_col==m.nblk_row);
  const int nblk = m.nblk_col;
  assert(m.len_col==m.len_row);
  const int len = m.len_col;

  assert(m.is_dia);
  mat.Initialize(nblk, len, true);

  m_diaInd.resize(nblk);

  for (int iblk = 0; iblk<nblk; ++iblk){
    std::vector<CRowLevNext> listNonzero;
    {	// copy row pattern of input matrix into listNonzero
      listNonzero.resize(m.colInd[iblk+1]-m.colInd[iblk]);
      int inz = 0;
      for (unsigned int ijcrs = m.colInd[iblk]; ijcrs<m.colInd[iblk+1]; ijcrs++){
        assert(ijcrs<m.ncrs);
        const int jblk0 = m.rowPtr[ijcrs];
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
      const int kblk0 = listNonzero[knz_cur].row;
      assert(kblk0<nblk);
      const int ik_lev0 = listNonzero[knz_cur].lev;
      if (ik_lev0+1>lev_fill && lev_fill!=-1){
        knz_cur = listNonzero[knz_cur].next;
        if (knz_cur==-1) break;
        continue;
      }
      if (kblk0>=iblk) break;

      int jnz_cur = knz_cur;
      for (unsigned int kjcrs = m_diaInd[kblk0]; kjcrs<mat.colInd[kblk0+1]; kjcrs++){
        const int kj_lev0 = aRowLev[kjcrs].lev;
        if (kj_lev0+1>lev_fill && lev_fill!=-1) continue;
        const int jblk0 = aRowLev[kjcrs].row;
        assert(jblk0>kblk0 && jblk0<nblk);
        assert(listNonzero[jnz_cur].row < jblk0);
        if (jblk0==iblk) continue; // already filled-in on the diagonal

        // check if this is fill in
        bool is_fill_in = false;
        for (;;){
          const int jnz_nex = listNonzero[jnz_cur].next;
          assert((jnz_nex>=0&&jnz_nex<nblk)||jnz_nex==-1);
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
      assert((knz_cur>=0&&knz_cur<nblk)||knz_cur==-1);
      if (knz_cur==-1) break;
    }

    ////////////////////

    {     
//      if (ColInd_pre[iblk]+nonzero.size() > m_pRowLev->max_size()){
//        std::cout<<"		overflow and memory reallocate in ilu frac"<<std::endl;
//      }
      // copy "listNonzero" to "aRowLev"
      aRowLev.resize(mat.colInd[iblk]+listNonzero.size());
      int icrs0 = mat.colInd[iblk];
      for (int inz = 0; inz!=-1; inz = listNonzero[inz].next){
        const int jblk = listNonzero[inz].row;
        const int jlev = listNonzero[inz].lev;
        assert(jblk<nblk);
        assert(jblk!=iblk);
        aRowLev[icrs0].row = jblk;
        aRowLev[icrs0].lev = jlev;
        icrs0++;
      }

      mat.colInd[iblk+1] = icrs0;
      mat.ncrs = icrs0;
      m_diaInd[iblk] = icrs0;
      for (unsigned int ijcrs = mat.colInd[iblk]; ijcrs<mat.colInd[iblk+1]; ijcrs++){
        const int jblk0 = aRowLev[ijcrs].row;
        if (jblk0 > iblk){
          m_diaInd[iblk] = ijcrs;
          break;
        }
      }
    }
  }

  {
    const unsigned int ncrs = mat.ncrs;
    assert(aRowLev.size()==ncrs);
    mat.rowPtr.resize(ncrs);
    for (unsigned int icrs = 0; icrs<ncrs; ++icrs){
      mat.rowPtr[icrs] = aRowLev[icrs].row;
    }
    /*
    for (int iblk = 0; iblk<nblk; ++iblk){
      for (int icrs = mat.m_colInd[iblk]; icrs<mat.m_colInd[iblk+1]; ++icrs){
        const int jblk = mat.m_rowPtr[icrs];
//        std::cout<<iblk<<" "<<jblk<<" "<<aRowLev[icrs].lev<<std::endl;
      }
    }
     */
    const int blksize = len*len;
    mat.valCrs = new double[ncrs*blksize];
    assert(mat.is_dia);
    mat.valDia = new double[nblk*blksize];
    for (int i = 0; i<nblk*blksize; i++){ mat.valDia[i] = m.valDia[i]; }
    std::cout<<"ncrs: "<<ncrs<<" "<<m.ncrs<<std::endl;
  }

}

void CPreconditionerILU::ForwardSubstitution( std::vector<double>& vec ) const
{
  const int len = mat.len_col;
  const int nblk = mat.nblk_col;
  
	if( len == 1 ){
		const unsigned int* colind = mat.colInd.data();
		const unsigned int* rowptr = mat.rowPtr.data();
		const double* vcrs = mat.valCrs;
		const double* vdia = mat.valDia;
		////////////////    
		for(int iblk=0;iblk<nblk;iblk++){	// ëOêiè¡ãé
			double lvec_i = vec[iblk];
			for(unsigned int ijcrs=colind[iblk];ijcrs<m_diaInd[iblk];ijcrs++){
				assert( ijcrs<mat.ncrs );
				const int jblk0 = rowptr[ijcrs];
				assert( jblk0<iblk );
				lvec_i -= vcrs[ijcrs]*vec[jblk0];
			}
			vec[iblk] = vdia[iblk]*lvec_i;
		}
	}
	else if( len == 2 ){
		const unsigned int* colind = mat.colInd.data();
		const unsigned int* rowptr = mat.rowPtr.data();
		const double* vcrs = mat.valCrs;
		const double* vdia = mat.valDia;
		////////////////
		double pTmpVec[2];
		for(int iblk=0;iblk<nblk;iblk++){
			pTmpVec[0] = vec[iblk*2+0];
			pTmpVec[1] = vec[iblk*2+1];
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = m_diaInd[iblk];
			for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
				assert( ijcrs<mat.ncrs );
				const int jblk0 = rowptr[ijcrs];
				assert( jblk0<iblk );
				const double* vij = &vcrs[ijcrs*4];
        const double valj0 = vec[jblk0*2+0];
				const double valj1 = vec[jblk0*2+1];
				pTmpVec[0] -= vij[0]*valj0+vij[1]*valj1;
				pTmpVec[1] -= vij[2]*valj0+vij[3]*valj1;
			}
			const double* vii = &vdia[iblk*4];
			vec[iblk*2+0] = vii[0]*pTmpVec[0]+vii[1]*pTmpVec[1];
			vec[iblk*2+1] = vii[2]*pTmpVec[0]+vii[3]*pTmpVec[1];
		}
	}
	else if( len == 3 ){
		const unsigned int* colind = mat.colInd.data();
		const unsigned int* rowptr = mat.rowPtr.data();
		const double* vcrs = mat.valCrs;
		const double* vdia = mat.valDia;
		////////////////
		double pTmpVec[3];
		for(int iblk=0;iblk<nblk;iblk++){
			pTmpVec[0] = vec[iblk*3+0];
			pTmpVec[1] = vec[iblk*3+1];
			pTmpVec[2] = vec[iblk*3+2];
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = m_diaInd[iblk];
			for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
				assert( ijcrs<mat.ncrs );
				const int jblk0 = rowptr[ijcrs];
				assert( jblk0<iblk );
				const double* vij = &vcrs[ijcrs*9];
				const double valj0 = vec[jblk0*3+0];
				const double valj1 = vec[jblk0*3+1];
				const double valj2 = vec[jblk0*3+2];
				pTmpVec[0] -= vij[0]*valj0+vij[1]*valj1+vij[2]*valj2;
				pTmpVec[1] -= vij[3]*valj0+vij[4]*valj1+vij[5]*valj2;
				pTmpVec[2] -= vij[6]*valj0+vij[7]*valj1+vij[8]*valj2;
			}
			const double* vii = &vdia[iblk*9];
			vec[iblk*3+0] = vii[0]*pTmpVec[0]+vii[1]*pTmpVec[1]+vii[2]*pTmpVec[2];
			vec[iblk*3+1] = vii[3]*pTmpVec[0]+vii[4]*pTmpVec[1]+vii[5]*pTmpVec[2];
			vec[iblk*3+2] = vii[6]*pTmpVec[0]+vii[7]*pTmpVec[1]+vii[8]*pTmpVec[2];
		}
	}
  else if (len==4){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const double* vcrs = mat.valCrs;
    const double* vdia = mat.valDia;
    ////////////////
    double pTmpVec[4];
    for (int iblk = 0; iblk<nblk; iblk++){
      pTmpVec[0] = vec[iblk*4+0];
      pTmpVec[1] = vec[iblk*4+1];
      pTmpVec[2] = vec[iblk*4+2];
      pTmpVec[3] = vec[iblk*4+3];
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = m_diaInd[iblk];
      for (unsigned int ijcrs = icrs0; ijcrs<icrs1; ijcrs++){
        assert(ijcrs<mat.ncrs);
        const int jblk0 = rowptr[ijcrs];
        assert(jblk0<iblk);
        const double* vij = &vcrs[ijcrs*16];
        const double valj0 = vec[jblk0*4+0];
        const double valj1 = vec[jblk0*4+1];
        const double valj2 = vec[jblk0*4+2];
        const double valj3 = vec[jblk0*4+3];
        pTmpVec[0] -= vij[ 0]*valj0+vij[ 1]*valj1+vij[ 2]*valj2+vij[ 3]*valj3;
        pTmpVec[1] -= vij[ 4]*valj0+vij[ 5]*valj1+vij[ 6]*valj2+vij[ 7]*valj3;
        pTmpVec[2] -= vij[ 8]*valj0+vij[ 9]*valj1+vij[10]*valj2+vij[11]*valj3;
        pTmpVec[3] -= vij[12]*valj0+vij[13]*valj1+vij[14]*valj2+vij[15]*valj3;
      }
      const double* vii = &vdia[iblk*16];
      vec[iblk*4+0] = vii[ 0]*pTmpVec[0]+vii[ 1]*pTmpVec[1]+vii[ 2]*pTmpVec[2]+vii[ 3]*pTmpVec[3];
      vec[iblk*4+1] = vii[ 4]*pTmpVec[0]+vii[ 5]*pTmpVec[1]+vii[ 6]*pTmpVec[2]+vii[ 7]*pTmpVec[3];
      vec[iblk*4+2] = vii[ 8]*pTmpVec[0]+vii[ 9]*pTmpVec[1]+vii[10]*pTmpVec[2]+vii[11]*pTmpVec[3];
      vec[iblk*4+3] = vii[12]*pTmpVec[0]+vii[13]*pTmpVec[1]+vii[14]*pTmpVec[2]+vii[15]*pTmpVec[3];
    }
  }
	else{
		const int blksize = len*len;
		double* pTmpVec = new double [len];
		for(int iblk=0;iblk<nblk;iblk++){
			for(int idof=0;idof<len;idof++){
				pTmpVec[idof] = vec[iblk*len+idof];
			}
			for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<m_diaInd[iblk];ijcrs++){
				assert( ijcrs<mat.ncrs );
				const int jblk0 = mat.rowPtr[ijcrs];
				assert( jblk0<iblk );
				const double* vij = &mat.valCrs[ijcrs*blksize];
				for(int idof=0;idof<len;idof++){
					for(int jdof=0;jdof<len;jdof++){
						pTmpVec[idof] -= vij[idof*len+jdof]*vec[jblk0*len+jdof];
					}
				}
			}
			const double* vii = &mat.valDia[iblk*blksize];
			for(int idof=0;idof<len;idof++){
				double dtmp1 = 0.0;
				for(int jdof=0;jdof<len;jdof++){
					dtmp1 += vii[idof*len+jdof]*pTmpVec[jdof];
				}
				vec[iblk*len+idof] = dtmp1;
			}
		}
		delete[] pTmpVec;
	}
}

void CPreconditionerILU::BackwardSubstitution( std::vector<double>& vec ) const
{
  const int len = mat.len_col;
  const int nblk = mat.nblk_col;

	if( len == 1 ){
		const unsigned int* colind = mat.colInd.data();
		const unsigned int* rowptr = mat.rowPtr.data();
		const double* vcrs = mat.valCrs;
		////////////////
		for(int iblk=nblk-1;iblk>=0;iblk--){	
			assert( (int)iblk < nblk );
			double lvec_i = vec[iblk];
			for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<colind[iblk+1];ijcrs++){
				assert( ijcrs<mat.ncrs );
				const int jblk0 = rowptr[ijcrs];
//        std::cout << jblk0 << " " << iblk << " " << nblk << std::endl;
				assert( jblk0>(int)iblk && jblk0<nblk );
				lvec_i -=  vcrs[ijcrs]*vec[jblk0];
			}
			vec[iblk] = lvec_i;
		}
	}
	else if( len == 2 ){
		const unsigned int* colind = mat.colInd.data();
		const unsigned int* rowptr = mat.rowPtr.data();
		const double* vcrs = mat.valCrs;
    ////////////////
		double pTmpVec[2];
		for(int iblk=nblk-1;iblk>=0;iblk--){
			assert( (int)iblk < nblk );
			pTmpVec[0] = vec[iblk*2+0];
			pTmpVec[1] = vec[iblk*2+1];
			const unsigned int icrs0 = m_diaInd[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
				assert( ijcrs<mat.ncrs );
				const int jblk0 = rowptr[ijcrs];
				assert( jblk0>(int)iblk && jblk0<nblk );
				const double* vij = &vcrs[ijcrs*4];
				const double valj0 = vec[jblk0*2+0];
				const double valj1 = vec[jblk0*2+1];
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
		const double* vcrs = mat.valCrs;
    ////////////////
		double pTmpVec[3];
		for(int iblk=nblk-1;iblk>=0;iblk--){
			assert( (int)iblk < nblk );
			pTmpVec[0] = vec[iblk*3+0];
			pTmpVec[1] = vec[iblk*3+1];
			pTmpVec[2] = vec[iblk*3+2];
			const int icrs0 = m_diaInd[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
				assert( ijcrs<mat.ncrs );
				const int jblk0 = rowptr[ijcrs];
				assert( jblk0>(int)iblk && jblk0<nblk );
				const double* vij = &vcrs[ijcrs*9];
				const double valj0 = vec[jblk0*3+0];
				const double valj1 = vec[jblk0*3+1];
				const double valj2 = vec[jblk0*3+2];
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
    const double* vcrs = mat.valCrs;
    ////////////////
    double pTmpVec[4];
    for (int iblk = nblk-1; iblk>=0; iblk--){
      assert((int)iblk < nblk);
      pTmpVec[0] = vec[iblk*4+0];
      pTmpVec[1] = vec[iblk*4+1];
      pTmpVec[2] = vec[iblk*4+2];
      pTmpVec[3] = vec[iblk*4+3];
      const int icrs0 = m_diaInd[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for (unsigned int ijcrs = icrs0; ijcrs<icrs1; ijcrs++){
        assert(ijcrs<mat.ncrs);
        const int jblk0 = rowptr[ijcrs];
        assert(jblk0>(int)iblk && jblk0<nblk);
        const double* vij = &vcrs[ijcrs*16];
        const double valj0 = vec[jblk0*4+0];
        const double valj1 = vec[jblk0*4+1];
        const double valj2 = vec[jblk0*4+2];
        const double valj3 = vec[jblk0*4+3];
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
		double* pTmpVec = new double [len];	// çÏã∆ópÇÃè¨Ç≥Ç»îzóÒ
		for(int iblk=nblk-1;iblk>=0;iblk--){
			assert( (int)iblk < nblk );
			for(int idof=0;idof<len;idof++){
				pTmpVec[idof] = vec[iblk*len+idof];
			}
			for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
				assert( ijcrs<mat.ncrs );
				const int jblk0 = mat.rowPtr[ijcrs];
				assert( jblk0>(int)iblk && jblk0<nblk );
				const double* vij = &mat.valCrs[ijcrs*blksize];
				for(int idof=0;idof<len;idof++){
					for(int jdof=0;jdof<len;jdof++){
						pTmpVec[idof] -= vij[idof*len+jdof]*vec[jblk0*len+jdof];
					}
				}
			}
			for(int idof=0;idof<len;idof++){
				vec[iblk*len+idof] = pTmpVec[idof];
			}
		}
		delete[] pTmpVec;
	}
}

void CPreconditionerILU::SetValueILU(const CMatrixSparse& m)
{
  const int nblk = mat.nblk_col;
	const int len = mat.len_col;
	const int blksize = len*len;
//  for(int i=0;i<mat.m_ncrs*blksize;i++){ mat.m_valCrs[i] = m.m_valCrs[i]; }
  std::vector<int> row2crs(nblk,-1);
  for(int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<mat.ncrs );
      const int jblk0 = mat.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = ijcrs;
    }
    for(unsigned int ijcrs=m.colInd[iblk];ijcrs<m.colInd[iblk+1];ijcrs++){
      assert( ijcrs<m.ncrs );
      const int jblk0 = m.rowPtr[ijcrs];
      assert( jblk0<nblk );
      const int ijcrs0 = row2crs[jblk0];
      if( ijcrs0 == -1 ) continue;
      const double* pval_in = &m.valCrs[ijcrs*blksize];
      double* pval_out = &mat.valCrs[ijcrs0*blksize];
      for(int i=0;i<blksize;i++){ *(pval_out+i) = *(pval_in+i); }
    }
    for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<mat.ncrs );
      const int jblk0 = mat.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = -1;
    }
  }
  for(int i=0;i<nblk*blksize;i++){ mat.valDia[i] = m.valDia[i]; }
}

// numerical factorization
bool CPreconditionerILU::DoILUDecomp()
{
  const int nmax_sing = 10;
	int icnt_sing = 0;
  
	const int len = mat.len_col;
  const int nblk = mat.nblk_col;
//  const int m_ncrs = mat.m_ncrs;
  const unsigned int* colind = mat.colInd.data();
  const unsigned int* rowptr = mat.rowPtr.data();
  double* vcrs = mat.valCrs;
  double* vdia = mat.valDia;
  const unsigned int m_ncrs = colind[nblk];
  
  std::vector<int> row2crs(nblk,-1);
  
	if( len == 1 ){
		for(int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const int jblk0 = rowptr[ijcrs];
        assert( jblk0<nblk );
				row2crs[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U] 
			for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
				const int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double ikvalue = vcrs[ikcrs];
				for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
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
				const int jblk0 = rowptr[ijcrs];
        assert( jblk0<nblk );
				row2crs[jblk0] = -1;
			}
		}	// end iblk
	}
	////////////////////////////////////////////////////////////////
	else if( len == 2 ){
		double TmpBlk[4];
		for(int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U]
			for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
				const int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double* vik = &vcrs[ikcrs*4];
				for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
					double* vkj = &vcrs[kjcrs*4]; assert( vkj != 0 );
					double* vij = 0;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs[jblk0];
						if( ijcrs0 == -1 ) continue;	
						vij = &vcrs[ijcrs0*4];
					}
          else{ vij = &vdia[iblk*4]; }
          assert( vij != 0 );
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
				const int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = -1;
			}
		}	// end iblk
	}
	////////////////////////////////////////////////////////////////
	else if( len == 3 ){	// lenBlk >= 3
		double tmpBlk[9];
		for(int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U]
			for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
				const int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double* vik = &vcrs[ikcrs*9];
				for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
					double* vkj = &vcrs[kjcrs*9]; assert( vkj != 0 );
					double* vij = 0;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs[jblk0];
            if( ijcrs0 == -1 ){ continue; }
						vij = &vcrs[ijcrs0*9];
					}
          else{ vij = &vdia[iblk*9]; }
					assert( vij != 0 );
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
				const int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = -1;
			}
		}	// end iblk
	}
  ////////////////////////////////////////////////////////////////
	else{	// lenBlk >= 4
    const int blksize = len*len;
		double* pTmpBlk = new double [blksize];
		for(int iblk=0;iblk<nblk;iblk++){
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){
        assert( ijcrs<m_ncrs );
				const int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U] 
			for(unsigned int ikcrs=colind[iblk];ikcrs<m_diaInd[iblk];ikcrs++){
				const int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double* vik = &vcrs[ikcrs*blksize];
				for(unsigned int kjcrs=m_diaInd[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
					double* vkj = &vcrs[kjcrs *blksize]; assert( vkj != 0 );
					double* vij = 0;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs[jblk0];
            if( ijcrs0 == -1 ){ continue; }
						vij = &vcrs[ijcrs0*blksize];
					}
          else{ vij = &vdia[iblk *blksize]; }
					assert( vij != 0 );
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
				const int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs[jblk0] = -1;
			}
		}	// end iblk
		delete[] pTmpBlk;
	}  
	return true;
}



std::vector<double> Solve_PCG
(double* r_vec,
 double* x_vec,
 double conv_ratio,
 int iteration,
 const CMatrixSparse& mat,
 const CPreconditionerILU& ilu)
{
	const int nblk = mat.nblk_col;
  const int len = mat.len_col;
//  assert(r_vec.size() == nblk*len);
  const int ndof = nblk*len;
  std::vector<double> aResHistry;
  
  // {x} = 0
//  x_vec.resize(ndof);
  for(int i=0;i<ndof;i++){ x_vec[i] = 0; }
  
	double inv_sqnorm_res0;
	{
		const double sqnorm_res0 = InnerProduct(r_vec,r_vec,ndof);
    aResHistry.push_back(sqnorm_res0);
		if( sqnorm_res0 < 1.0e-30 ){ return aResHistry; }
		inv_sqnorm_res0 = 1.0 / sqnorm_res0;
	}
  
  // {Pr} = [P]{r}
  std::vector<double> Pr_vec(r_vec,r_vec+ndof);
  ilu.Solve(Pr_vec);
  // {p} = {Pr}
  std::vector<double> p_vec = Pr_vec;
  
  // rPr = ({r},{Pr})
	double rPr = InnerProduct(r_vec,Pr_vec.data(),ndof);
	for(int iitr=0;iitr<iteration;iitr++){
    
		{
      std::vector<double>& Ap_vec = Pr_vec;      
      // {Ap} = [A]{p}
			mat.MatVec(1.0,p_vec,0.0,Ap_vec);
      // alpha = ({r},{Pr})/({p},{Ap})
			const double pAp = InnerProduct(p_vec,Ap_vec);
			double alpha = rPr / pAp;
      // {r} = -alpha*{Ap} + {r}
      AXPY(-alpha,Ap_vec.data(),r_vec, ndof);
      // {x} = +alpha*{p } + {x}
      AXPY(+alpha,p_vec.data(), x_vec, ndof);
    }
    
		{	// Converge Judgement
			double sqnorm_res = InnerProduct(r_vec,r_vec,ndof);
      aResHistry.push_back(sqrt(sqnorm_res));
			if( sqnorm_res * inv_sqnorm_res0 < conv_ratio*conv_ratio ){
				return aResHistry;
			}
		}
    
		{	// calc beta
      // {Pr} = [P]{r}
      for(int i=0;i<ndof;i++){ Pr_vec[i] = r_vec[i]; }
			ilu.Solve(Pr_vec);
      // rPr1 = ({r},{Pr})
			const double rPr1 = InnerProduct(r_vec,Pr_vec.data(),ndof);
      // beta = rPr1/rPr
			double beta = rPr1/rPr;
			rPr = rPr1;
      // {p} = {Pr} + beta*{p}
      for(int i=0;i<ndof;i++){ p_vec[i] = Pr_vec[i] + beta*p_vec[i]; }
    }
	}
  {
    // Converge Judgement
    double sq_norm_res = InnerProduct(r_vec,r_vec,ndof);
    aResHistry.push_back(sqrt(sq_norm_res));
  }
  return aResHistry;
}


std::vector<double> Solve_PBiCGStab
(double* r_vec,
 double* x_vec,
 double conv_ratio,
 int num_iter,
 const CMatrixSparse& mat,
 const CPreconditionerILU& ilu)
{
  assert( mat.is_dia );
  assert( mat.nblk_col == mat.nblk_row );
  assert( mat.len_col == mat.len_row );
  
  const int nblk = mat.nblk_col;
  const int len = mat.len_col;
//  assert(r_vec.size() == nblk*len);
  const int ndof = nblk*len;
  std::vector<double> aResHistry;
  
  const unsigned int max_iter = num_iter;
  const double tolerance = conv_ratio;
 
  // {u} = 0
  for(int i=0;i<ndof;++i){ x_vec[i] = 0.0; }
  
  double sq_inv_norm_res_ini;
  {
    const double sq_norm_res_ini = InnerProduct(r_vec,r_vec,ndof);
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
  std::vector<double> p_vec(ndof);
  std::vector<double> Mp_vec(ndof);
  std::vector<double> AMp_vec(ndof);
  std::vector<double> r2_vec(ndof);
  
  // {r2} = {r}
  r2_vec.assign(r_vec,r_vec+ndof);
  
  // {p} = {r}
  p_vec.assign(r_vec,r_vec+ndof);
  
  num_iter = max_iter;
  for(unsigned int iitr=1;iitr<max_iter;iitr++)
  {
    // {Mp_vec} = [M^-1]*{p}
//    ls.COPY(ip,iMp);
    Mp_vec = p_vec;
//    ls.SolvePrecond(iMp);
    ilu.Solve(Mp_vec);
    
    // calc (r,r0*)
    const double r_r2 = InnerProduct(r_vec,r2_vec.data(),ndof);
    
    //        std::cout << "r_r2 : " << r_r2 << std::endl;
    
    // calc {AMp_vec} = [A]*{Mp_vec}
    mat.MatVec(1.0, Mp_vec, 0.0, AMp_vec);
    
    // calc alpha
    double alpha;
    {
      const double denominator = InnerProduct(AMp_vec,r2_vec);
      alpha = r_r2 / denominator;
    }
    
    //        std::cout << "Alpha : " << alpha << std::endl;
    
    // calc s_vector
    s_vec.assign(r_vec,r_vec+ndof);
    AXPY(-alpha,AMp_vec,s_vec);
//    ls.COPY(ir,is);
//    ls.AXPY(-alpha,iAMp,is);
    
    //        std::cout << "ir iAMp is " << ls.DOT(ir,ir) << " " << ls.DOT(iAMp,iAMp) << " " << ls.DOT(is,is) << std::endl;
    
    // {Ms_vec} = [M^-1]*{s}
//    ls.COPY(is,iMs);
    Ms_vec = s_vec;
//    ls.SolvePrecond(iMs);
    ilu.Solve(Ms_vec);
    
    //        std::cout << "Is iMs " << ls.DOT(is,is) << " " << ls.DOT(iMs,iMs) << std::endl;
    
    // calc {AMs_vec} = [A]*{Ms_vec}
//    ls.MATVEC(1.0,iMs,0.0,iAMs);
    mat.MatVec(1.0,Ms_vec,0.0,AMs_vec);
    
    double omega;
    {	// calc omega
      const double denominator = InnerProduct(AMs_vec,AMs_vec);
      const double numerator = InnerProduct(s_vec,AMs_vec);
      //            std::cout << "Omega0 : " << denominator << " " << numerator << std::endl;
      omega = numerator / denominator;
    }
    
    //        std::cout << "Omega : " << omega << std::endl;
    
    // update solution
//    ls.AXPY(alpha,iMp,ix);
    AXPY(alpha,Mp_vec.data(),x_vec,ndof);
//    ls.AXPY(omega,iMs,ix);
    AXPY(omega,Ms_vec.data(),x_vec,ndof);
    
    // update residual
//    ls.COPY(is,ir);
    for(int i=0;i<ndof;++i){ r_vec[i] = s_vec[i]; }
    
//    ls.AXPY(-omega,iAMs,ir);
    AXPY(-omega,AMs_vec.data(),r_vec,ndof);
    
    {
      const double sq_norm_res = InnerProduct(r_vec,r_vec,ndof);
      const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res_ini;
//      std::cout << iitr << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
      if( sq_conv_ratio < tolerance * tolerance ){
        aResHistry.push_back( sqrt( sq_norm_res ) );
        return aResHistry;
      }
    }
    
    double beta;
    {	// calc beta
      const double tmp1 = InnerProduct(r_vec,r2_vec.data(),ndof);
      beta = tmp1 * alpha / (r_r2*omega);
    }
    
    // update p_vector
//    ls.SCAL(beta,ip);
    for(int i=0;i<ndof;++i){ p_vec[i] *= beta; }
//    ls.AXPY(1.0,ir,ip);
    AXPY(1.0,r_vec,p_vec.data(),ndof);
//    ls.AXPY(-beta*omega,iAMp,ip);
    AXPY(-beta*omega,AMp_vec,p_vec);
  }
  
  return aResHistry;
}

