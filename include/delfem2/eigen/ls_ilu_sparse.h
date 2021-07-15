#ifndef DFM2_EIGEN_LS_ILU_SPARSE_H
#define DFM2_EIGEN_LS_ILU_SPARSE_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <climits>
#include "delfem2/eigen/ls_sparse.h"

namespace delfem2 {

template<class MAT, class ALLOCATOR>
class CILU_SparseBlock {
public:
  CILU_SparseBlock() {}
public:
  unsigned int nblk;
  std::vector<unsigned int> colInd;
  std::vector<unsigned int> rowPtr;
  std::vector<MAT, ALLOCATOR> valCrs;
  std::vector<MAT, ALLOCATOR> valDia;
  std::vector<unsigned int> m_diaInd;
};

// -------------------------

template<class MAT, class ALLOCATOR>
void ILU_SetPattern0(
    CILU_SparseBlock<MAT,ALLOCATOR>& ilu,
    const CMatrixSparseBlock<MAT,ALLOCATOR>& A)
{
  const unsigned int nblk = A.nrowblk;
  ilu.nblk = nblk;
  ilu.colInd = A.colInd;
  ilu.rowPtr = A.rowPtr;
  ilu.valCrs.resize(ilu.rowPtr.size());
  ilu.valDia.resize(ilu.nblk);
  ilu.m_diaInd.resize(nblk);
  // ---------------
  // sort mat.rowPtr
  for(unsigned int iblk=0;iblk<nblk;++iblk) {
    const unsigned int icrs0 = ilu.colInd[iblk];
    const unsigned int icrs1 = ilu.colInd[iblk+1];
    std::sort(
        ilu.rowPtr.data()+icrs0,
        ilu.rowPtr.data()+icrs1);
  }
  // ------------
  // set m_diaInd
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    ilu.m_diaInd[iblk] = ilu.colInd[iblk+1];
    for(unsigned int icrs=ilu.colInd[iblk];icrs<ilu.colInd[iblk+1];icrs++){
      assert( icrs < ilu.rowPtr.size() );
      const int jblk0 = ilu.rowPtr[icrs];
      assert( jblk0 < nblk );
      if( jblk0 > iblk ){
        ilu.m_diaInd[iblk] = icrs;
        break;
      }
    }
  }
}

template<class MAT, class ALLOCATOR>
void ILU_CopyValue(
    CILU_SparseBlock<MAT,ALLOCATOR>& ilu,
    const CMatrixSparseBlock<MAT,ALLOCATOR>& A)
{
  const unsigned int nblk = ilu.nblk;
  assert( A.nrowblk == nblk && A.ncolblk == nblk );
  std::vector<int> row2crs(nblk,-1);
  // copy diagonal value
  ilu.valDia = A.valDia;
  // copy off-diagonal values
  for(auto& b : ilu.valCrs){ b.setZero(); } // set zero
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=ilu.colInd[iblk];ijcrs<ilu.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      const unsigned int jblk0 = ilu.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = ijcrs;
    }
    for(unsigned int ijcrs=A.colInd[iblk];ijcrs<A.colInd[iblk+1];ijcrs++){
      assert( ijcrs<A.rowPtr.size() );
      const unsigned int jblk0 = A.rowPtr[ijcrs];
      assert( jblk0<nblk );
      const int ijcrs0 = row2crs[jblk0];
      if( ijcrs0 == UINT_MAX ) continue;
      ilu.valCrs[ijcrs0] = A.valCrs[ijcrs];
    }
    for(unsigned int ijcrs=ilu.colInd[iblk];ijcrs<ilu.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      const unsigned int jblk0 = ilu.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = -1;
    }
  }
}

template<class MAT, class ALLOCATOR>
bool ILU_Decompose(
    CILU_SparseBlock<MAT,ALLOCATOR>& ilu)
{
  const unsigned int nblk = ilu.nblk;
  std::vector<unsigned int> row2crs(nblk,UINT_MAX);
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=ilu.colInd[iblk];ijcrs<ilu.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      const unsigned int jblk0 = ilu.rowPtr[ijcrs];
      assert( jblk0<nblk );
      row2crs[jblk0] = ijcrs;
    }
    // [L] * [D^-1*U]
    for(unsigned int ikcrs=ilu.colInd[iblk];ikcrs<ilu.m_diaInd[iblk];ikcrs++){
      const unsigned int kblk = ilu.rowPtr[ikcrs];
      assert( kblk<nblk );
      const MAT& vik = ilu.valCrs[ikcrs];
      for(unsigned int kjcrs=ilu.m_diaInd[kblk];kjcrs<ilu.colInd[kblk+1];kjcrs++){
        const unsigned int jblk0 = ilu.rowPtr[kjcrs];
        assert( jblk0<nblk );
        const MAT& vkj = ilu.valCrs[kjcrs];
        if( jblk0 != iblk ){
          const unsigned int ijcrs0 = row2crs[jblk0];
          if( ijcrs0 == UINT_MAX ){ continue; }
//          ilu.mat.valCrs[ijcrs0] -= vik*vkj.transpose();
//          ilu.mat.valCrs[ijcrs0] -= vik.transpose()*vkj;
          ilu.valCrs[ijcrs0] -= vik*vkj;
        }
        else{
//          ilu.mat.valDia[iblk] -= vik*vkj.transpose();
//          ilu.mat.valDia[iblk] -= vik.transpose()*vkj;
          ilu.valDia[iblk] -= vik*vkj;
        }
      }
    }
    // invserse diagonal
    ilu.valDia[iblk] = ilu.valDia[iblk].inverse().eval();
    // [U] = [1/D][U]
    for(unsigned int ijcrs=ilu.m_diaInd[iblk];ijcrs<ilu.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      ilu.valCrs[ijcrs] = ilu.valDia[iblk]*ilu.valCrs[ijcrs];
    }
    for(unsigned int ijcrs=ilu.colInd[iblk];ijcrs<ilu.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      const unsigned int jblk0 = ilu.rowPtr[ijcrs];
      assert( jblk0<nblk );
      row2crs[jblk0] = UINT_MAX;
    }
  }	// end iblk
  return true;
}

template<typename REAL, class MAT, class ALLOCATOR>
void SolvePrecond(
    Eigen::Matrix<REAL,-1,1>& vec,
    const CILU_SparseBlock<MAT,ALLOCATOR>& ilu)
{
  constexpr unsigned int nrowdim = MAT::RowsAtCompileTime;
  constexpr unsigned int ncoldim = MAT::ColsAtCompileTime;
  static_assert(nrowdim == ncoldim,
      "The block matrix need to be square");
  constexpr unsigned int ndim = nrowdim;
  // --------
  // forward
  const unsigned int nblk = ilu.nblk;
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=ilu.colInd[iblk];ijcrs<ilu.m_diaInd[iblk];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      const unsigned int jblk0 = ilu.rowPtr[ijcrs];
      assert( (int)jblk0<iblk );
      vec.template segment<ndim>(iblk*ndim) -= ilu.valCrs[ijcrs]*vec.template segment<ndim>(jblk0*ndim); // jblk0!=iblk
    }
    vec.template segment<ndim>(iblk*ndim) = ilu.valDia[iblk]*vec.template segment<ndim>(iblk*ndim).eval();
  }
  // -----
  // backward
  for(int iblk=nblk-1;iblk>=0;iblk--){
    assert( iblk < (int)nblk );
    for(unsigned int ijcrs=ilu.m_diaInd[iblk];ijcrs<ilu.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      const unsigned int jblk0 = ilu.rowPtr[ijcrs];
      assert( (int)jblk0>iblk && jblk0<nblk );
      vec.template segment<ndim>(iblk*ndim) -= ilu.valCrs[ijcrs]*vec.template segment<ndim>(jblk0*ndim); // jblk0!=iblk
    }
  }
}

template<typename REAL, class MAT, class ALLOCATOR>
void SolvePrecond(
    Eigen::Matrix<REAL,-1,MAT::RowsAtCompileTime,Eigen::RowMajor,-1,MAT::RowsAtCompileTime>& vec,
    const CILU_SparseBlock<MAT,ALLOCATOR>& ilu)
{
  constexpr unsigned int nrowdim = MAT::RowsAtCompileTime;
  constexpr unsigned int ncoldim = MAT::ColsAtCompileTime;
  static_assert(nrowdim == ncoldim,
      "The block matrix need to be square");
  // --------
  // forward
  const unsigned int nblk = ilu.nblk;
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=ilu.colInd[iblk];ijcrs<ilu.m_diaInd[iblk];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      const unsigned int jblk0 = ilu.rowPtr[ijcrs];
      assert( (int)jblk0<iblk );
      vec.row(iblk) -= ilu.valCrs[ijcrs]*vec.row(jblk0).transpose(); // jblk0!=iblk
    }
    vec.row(iblk) = ilu.valDia[iblk]*vec.row(iblk).transpose().eval();
  }
  // -----
  // backward
  for(int iblk=nblk-1;iblk>=0;iblk--){
    assert( iblk < (int)nblk );
    for(unsigned int ijcrs=ilu.m_diaInd[iblk];ijcrs<ilu.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.rowPtr.size() );
      const unsigned int jblk0 = ilu.rowPtr[ijcrs];
      assert( (int)jblk0>iblk && jblk0<nblk );
      vec.row(iblk) -= ilu.valCrs[ijcrs]*vec.row(jblk0).transpose(); // jblk0!=iblk
    }
  }
}

}

#endif