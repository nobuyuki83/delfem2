/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef ILU_SPARSE
#define ILU_SPARSE

#include <iostream>

#include "mats.h"

template <typename T>
class CPreconditionerILU
{
public:
  CPreconditionerILU(){}
  CPreconditionerILU(const CPreconditionerILU&); // copy
  ~CPreconditionerILU(){ m_diaInd.clear(); }
  void Initialize_ILU0(const CMatrixSparse<T>& m);
  void Initialize_ILUk(const CMatrixSparse<T>& m, int fill_level);
  void SetValueILU(const CMatrixSparse<T>& m);
  void Solve(std::vector<T>& vec) const{
		this->ForwardSubstitution(vec);
		this->BackwardSubstitution(vec);
  }
  bool DoILUDecomp();
private:
  void ForwardSubstitution(  std::vector<T>& vec ) const;
  void BackwardSubstitution( std::vector<T>& vec ) const;
public:
  CMatrixSparse<T> mat;
  std::vector<unsigned int> m_diaInd;
};




template <typename T>
std::vector<double> Solve_PCG(T* r_vec,
                              T* u_vec,
                              double conv_ratio,
                              int iteration,
                              const CMatrixSparse<T>& mat,
                              const CPreconditionerILU<T>& ilu);

template <typename T>
std::vector<double> Solve_PBiCGStab(T* r_vec,
                                    T* x_vec,
                                    double conv_ratio,
                                    int num_iter,
                                    const CMatrixSparse<T>& mat,
                                    const CPreconditionerILU<T>& ilu);

std::vector<double> Solve_PCOCG(std::complex<double>* r_vec,
                                std::complex<double>* x_vec,
                                double conv_ratio_tol,
                                int max_niter,
                                const CMatrixSparse<std::complex<double> >& mat,
                                const CPreconditionerILU<std::complex<double> >& ilu);

/*
template <typename T>
void SolveLinSys_PCG
(const CMatrixSparse<T>& mat_A,
 std::vector<double>& vec_b,
 std::vector<double>& vec_x,
 CPreconditionerILU<T>& ilu_A,
 double& conv_ratio,
 int& iteration)
{
  // set ILU preconditioner
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  // solve linear system
  //Solve_CG(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_BiCGSTAB(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_PBiCGSTAB(conv_ratio, iteration, mat_A, ilu_A, vec_b, vec_x);
  vec_x.resize(vec_b.size());
  Solve_PCG(vec_b.data(), vec_x.data(), conv_ratio, iteration,
            mat_A, ilu_A);
  std::cout<<"  conv_ratio:"<<conv_ratio<<"  iteration:"<<iteration<<std::endl;
}

template <typename T>
bool SolveLinSys_BiCGStab
(CMatrixSparse<T>& mat_A,
 std::vector<double>& vec_b,
 std::vector<double>& vec_x,
 CPreconditionerILU<T>& ilu_A,
 double& conv_ratio,
 int& iteration)
{
  // set ILU preconditioner
  ilu_A.SetValueILU(mat_A);
  bool res_ilu = ilu_A.DoILUDecomp();
  if( !res_ilu ){ return false; }
  // solve linear system
  //  double conv_ratio = 1.0e-4;
  //  int iteration = 1000;
  //  Solve_CG(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_BiCGSTAB(conv_ratio, iteration, mat_A, vec_b, vec_x);
  vec_x.resize(vec_b.size());
  Solve_PBiCGStab(vec_b.data(), vec_x.data(), conv_ratio, iteration, mat_A, ilu_A);
  /// Solve_PCG(conv_ratio, iteration, mat_A, ilu_A, vec_b, vec_x);
  //  std::cout<<"  interative solver --- conv_ratio:"<<conv_ratio<<"  iteration:"<<iteration<<std::endl;
  return true;
}
 */

template <typename T>
void CPreconditionerILU<T>::Initialize_ILU0
(const CMatrixSparse<T>& m)
{
  this->mat = m;
  const int nblk = m.nblk_col;
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

template <typename T>
void CPreconditionerILU<T>::SetValueILU
(const CMatrixSparse<T>& rhs)
{
  const int nblk = mat.nblk_col;
  const int len = mat.len_col;
  assert( rhs.nblk_col == nblk );
  assert( rhs.nblk_row == nblk );
  assert( rhs.len_col == len );
  assert( rhs.len_row == len );
  const int blksize = len*len;
  //  for(int i=0;i<mat.m_ncrs*blksize;i++){ mat.m_valCrs[i] = m.m_valCrs[i]; }
  std::vector<int> row2crs(nblk,-1);
  for(int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<mat.rowPtr.size() );
      const int jblk0 = mat.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = ijcrs;
    }
    for(unsigned int ijcrs=rhs.colInd[iblk];ijcrs<rhs.colInd[iblk+1];ijcrs++){
      assert( ijcrs<rhs.rowPtr.size() );
      const int jblk0 = rhs.rowPtr[ijcrs];
      assert( jblk0<nblk );
      const int ijcrs0 = row2crs[jblk0];
      if( ijcrs0 == -1 ) continue;
      const T* pval_in = &rhs.valCrs[ijcrs*blksize];
      T* pval_out = &mat.valCrs[ijcrs0*blksize];
      for(int i=0;i<blksize;i++){ *(pval_out+i) = *(pval_in+i); }
    }
    for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<mat.rowPtr.size() );
      const int jblk0 = mat.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = -1;
    }
  }
  for(int i=0;i<nblk*blksize;i++){ mat.valDia[i] = rhs.valDia[i]; }
}


template <typename T>
void CPreconditionerILU<T>::ForwardSubstitution
( std::vector<T>& vec ) const
{
  const int len = mat.len_col;
  const unsigned int nblk = mat.nblk_col;
  
  if( len == 1 ){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    const T* vdia = mat.valDia.data();
    ////////////////
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      T lvec_i = vec[iblk];
      for(unsigned int ijcrs=colind[iblk];ijcrs<m_diaInd[iblk];ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
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
    const T* vcrs = mat.valCrs.data();
    const T* vdia = mat.valDia.data();
    ////////////////
    T pTmpVec[2];
    for(int iblk=0;iblk<nblk;iblk++){
      pTmpVec[0] = vec[iblk*2+0];
      pTmpVec[1] = vec[iblk*2+1];
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = m_diaInd[iblk];
      for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = rowptr[ijcrs];
        assert( jblk0<iblk );
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
    ////////////////
    T pTmpVec[3];
    for(int iblk=0;iblk<nblk;iblk++){
      pTmpVec[0] = vec[iblk*3+0];
      pTmpVec[1] = vec[iblk*3+1];
      pTmpVec[2] = vec[iblk*3+2];
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = m_diaInd[iblk];
      for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = rowptr[ijcrs];
        assert( jblk0<iblk );
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
    ////////////////
    T pTmpVec[4];
    for (int iblk = 0; iblk<nblk; iblk++){
      pTmpVec[0] = vec[iblk*4+0];
      pTmpVec[1] = vec[iblk*4+1];
      pTmpVec[2] = vec[iblk*4+2];
      pTmpVec[3] = vec[iblk*4+3];
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = m_diaInd[iblk];
      for (unsigned int ijcrs = icrs0; ijcrs<icrs1; ijcrs++){
        assert(ijcrs<mat.rowPtr.size());
        const int jblk0 = rowptr[ijcrs];
        assert(jblk0<iblk);
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
    for(int iblk=0;iblk<nblk;iblk++){
      for(int idof=0;idof<len;idof++){
        pTmpVec[idof] = vec[iblk*len+idof];
      }
      for(unsigned int ijcrs=mat.colInd[iblk];ijcrs<m_diaInd[iblk];ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = mat.rowPtr[ijcrs];
        assert( jblk0<iblk );
        const T* vij = &mat.valCrs[ijcrs*blksize];
        for(int idof=0;idof<len;idof++){
          for(int jdof=0;jdof<len;jdof++){
            pTmpVec[idof] -= vij[idof*len+jdof]*vec[jblk0*len+jdof];
          }
        }
      }
      const T* vii = &mat.valDia[iblk*blksize];
      for(int idof=0;idof<len;idof++){
        T dtmp1 = 0.0;
        for(int jdof=0;jdof<len;jdof++){
          dtmp1 += vii[idof*len+jdof]*pTmpVec[jdof];
        }
        vec[iblk*len+idof] = dtmp1;
      }
    }
  }
}


template <typename T>
void CPreconditionerILU<T>::BackwardSubstitution
( std::vector<T>& vec ) const
{
  const unsigned int len = mat.len_col;
  const int nblk = mat.nblk_col;
  
  if( len == 1 ){
    const unsigned int* colind = mat.colInd.data();
    const unsigned int* rowptr = mat.rowPtr.data();
    const T* vcrs = mat.valCrs.data();
    ////////////////
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
    ////////////////
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
    ////////////////
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
    ////////////////
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
      for(int idof=0;idof<len;idof++){
        pTmpVec[idof] = vec[iblk*len+idof];
      }
      for(unsigned int ijcrs=m_diaInd[iblk];ijcrs<mat.colInd[iblk+1];ijcrs++){
        assert( ijcrs<mat.rowPtr.size() );
        const int jblk0 = mat.rowPtr[ijcrs];
        assert( jblk0>(int)iblk && jblk0<nblk );
        const T* vij = &mat.valCrs[ijcrs*blksize];
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
template <typename T>
void CPreconditionerILU<T>::Initialize_ILUk
(const CMatrixSparse<T>& m,
 int lev_fill)
{
  
  if (lev_fill==0){
    this->Initialize_ILU0(m);
    return;
  }
  
  std::vector<CRowLev> aRowLev;
  aRowLev.reserve(m.rowPtr.size()*4);
  
  assert(m.nblk_col==m.nblk_row);
  const int nblk = m.nblk_col;
  assert(m.len_col==m.len_row);
  const int len = m.len_col;
  
  assert(!m.valDia.empty());
  mat.Initialize(nblk, len, true);
  
  m_diaInd.resize(nblk);
  
  for(unsigned int iblk=0; iblk<nblk; ++iblk){
    std::vector<CRowLevNext> listNonzero;
    {  // copy row pattern of input matrix into listNonzero
      listNonzero.resize(m.colInd[iblk+1]-m.colInd[iblk]);
      int inz = 0;
      for (unsigned int ijcrs = m.colInd[iblk]; ijcrs<m.colInd[iblk+1]; ijcrs++){
        assert(ijcrs<m.rowPtr.size());
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
      mat.rowPtr.resize(icrs0);
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
    const unsigned int ncrs = mat.rowPtr.size();
    std::cout << aRowLev.size() << " " << ncrs << std::endl;
    assert(aRowLev.size()==ncrs);
    mat.rowPtr.resize(ncrs);
    for (unsigned int icrs = 0; icrs<ncrs; ++icrs){
      mat.rowPtr[icrs] = aRowLev[icrs].row;
    }
    const int blksize = len*len;
    mat.valCrs.resize(ncrs*blksize);
    assert(!mat.valDia.empty());
    mat.valDia = m.valDia;
    std::cout<<"ncrs: "<<ncrs<<" "<<m.rowPtr.size()<<std::endl;
  }
  
}

#endif 
