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
  CMatrixSparseBlock<MAT,ALLOCATOR> mat;
  std::vector<unsigned int> m_diaInd;
};

// -------------------------

template<class MAT, class ALLOCATOR>
void ILU_SetPattern0(
    CILU_SparseBlock<MAT,ALLOCATOR>& ilu,
    const CMatrixSparseBlock<MAT,ALLOCATOR>& A)
{
  ilu.mat.nrowblk = A.nrowblk;
  ilu.mat.ncolblk = A.ncolblk;
  ilu.mat.colInd = A.colInd;
  ilu.mat.rowPtr = A.rowPtr;
  ilu.mat.valCrs.resize(A.valCrs.size());
  ilu.mat.valDia.resize(A.valDia.size());
  const unsigned int nblk = A.nrowblk;
  ilu.m_diaInd.resize(nblk);
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    ilu.m_diaInd[iblk] = ilu.mat.colInd[iblk+1];
    for(unsigned int icrs=ilu.mat.colInd[iblk];icrs<ilu.mat.colInd[iblk+1];icrs++){
      assert( icrs < ilu.mat.rowPtr.size() );
      const int jblk0 = ilu.mat.rowPtr[icrs];
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
  const unsigned int nblk = ilu.mat.nrowblk;
  assert( A.nrowblk == nblk && A.ncolblk == nblk );
  std::vector<int> row2crs(nblk,-1);
  // copy diagonal value
  ilu.mat.valDia = A.valDia;
  // copy off-diagonal values
  for(auto& b : ilu.mat.valCrs){ b.setZero(); } // set zero
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=ilu.mat.colInd[iblk];ijcrs<ilu.mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.mat.rowPtr.size() );
      const unsigned int jblk0 = ilu.mat.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = ijcrs;
    }
    for(unsigned int ijcrs=A.colInd[iblk];ijcrs<A.colInd[iblk+1];ijcrs++){
      assert( ijcrs<A.rowPtr.size() );
      const unsigned int jblk0 = A.rowPtr[ijcrs];
      assert( jblk0<nblk );
      const int ijcrs0 = row2crs[jblk0];
      if( ijcrs0 == UINT_MAX ) continue;
      ilu.mat.valCrs[ijcrs0] = A.valCrs[ijcrs];
    }
    for(unsigned int ijcrs=ilu.mat.colInd[iblk];ijcrs<ilu.mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.mat.rowPtr.size() );
      const unsigned int jblk0 = ilu.mat.rowPtr[ijcrs];
      assert( jblk0 < nblk );
      row2crs[jblk0] = -1;
    }
  }
}

template<class MAT, class ALLOCATOR>
bool ILU_Decompose(
    CILU_SparseBlock<MAT,ALLOCATOR>& ilu)
{
  const unsigned int nblk = ilu.mat.nrowblk;
  std::vector<unsigned int> row2crs(nblk,UINT_MAX);
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=ilu.mat.colInd[iblk];ijcrs<ilu.mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.mat.rowPtr.size() );
      const unsigned int jblk0 = ilu.mat.rowPtr[ijcrs];
      assert( jblk0<nblk );
      row2crs[jblk0] = ijcrs;
    }
    // [L] * [D^-1*U]
    for(unsigned int ikcrs=ilu.mat.colInd[iblk];ikcrs<ilu.m_diaInd[iblk];ikcrs++){
      const unsigned int kblk = ilu.mat.rowPtr[ikcrs];
      assert( kblk<nblk );
      const MAT& vik = ilu.mat.valCrs[ikcrs];
      for(unsigned int kjcrs=ilu.m_diaInd[kblk];kjcrs<ilu.mat.colInd[kblk+1];kjcrs++){
        const unsigned int jblk0 = ilu.mat.rowPtr[kjcrs];
        assert( jblk0<nblk );
        const MAT& vkj = ilu.mat.valCrs[kjcrs];
        if( jblk0 != iblk ){
          const unsigned int ijcrs0 = row2crs[jblk0];
          if( ijcrs0 == UINT_MAX ){ continue; }
//          ilu.mat.valCrs[ijcrs0] -= vik*vkj.transpose();
//          ilu.mat.valCrs[ijcrs0] -= vik.transpose()*vkj;
          ilu.mat.valCrs[ijcrs0] -= vik*vkj;
        }
        else{
//          ilu.mat.valDia[iblk] -= vik*vkj.transpose();
//          ilu.mat.valDia[iblk] -= vik.transpose()*vkj;
          ilu.mat.valDia[iblk] -= vik*vkj;
        }
      }
    }
    // invserse diagonal
    ilu.mat.valDia[iblk] = ilu.mat.valDia[iblk].inverse().eval();
    // [U] = [1/D][U]
    for(unsigned int ijcrs=ilu.m_diaInd[iblk];ijcrs<ilu.mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.mat.rowPtr.size() );
      ilu.mat.valCrs[ijcrs] = ilu.mat.valDia[iblk]*ilu.mat.valCrs[ijcrs];
    }
    for(unsigned int ijcrs=ilu.mat.colInd[iblk];ijcrs<ilu.mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.mat.rowPtr.size() );
      const unsigned int jblk0 = ilu.mat.rowPtr[ijcrs];
      assert( jblk0<nblk );
      row2crs[jblk0] = UINT_MAX;
    }
  }	// end iblk
  return true;
}

template<typename REAL, class MAT, class ALLOCATOR>
void SolvePrecond(
    Eigen::VectorX<REAL>& vec,
    const CILU_SparseBlock<MAT,ALLOCATOR>& ilu)
{
  constexpr unsigned int nrowdim = MAT::RowsAtCompileTime;
  constexpr unsigned int ncoldim = MAT::ColsAtCompileTime;
  static_assert(nrowdim == ncoldim,
      "The block matrix need to be square");
  constexpr unsigned int ndim = nrowdim;
  // --------
  // forward
  const unsigned int nblk = ilu.mat.nrowblk;
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ijcrs=ilu.mat.colInd[iblk];ijcrs<ilu.m_diaInd[iblk];ijcrs++){
      assert( ijcrs<ilu.mat.rowPtr.size() );
      const unsigned int jblk0 = ilu.mat.rowPtr[ijcrs];
      assert( (int)jblk0<iblk );
      vec.template segment<ndim>(iblk*ndim) -= ilu.mat.valCrs[ijcrs]*vec.template segment<ndim>(jblk0*ndim); // jblk0!=iblk
    }
    vec.template segment<ndim>(iblk*ndim) = ilu.mat.valDia[iblk]*vec.template segment<ndim>(iblk*ndim).eval();
  }
  // -----
  // backward
  for(int iblk=nblk-1;iblk>=0;iblk--){
    assert( iblk < (int)nblk );
    for(unsigned int ijcrs=ilu.m_diaInd[iblk];ijcrs<ilu.mat.colInd[iblk+1];ijcrs++){
      assert( ijcrs<ilu.mat.rowPtr.size() );
      const unsigned int jblk0 = ilu.mat.rowPtr[ijcrs];
      assert( (int)jblk0>iblk && jblk0<nblk );
      vec.template segment<ndim>(iblk*ndim) -= ilu.mat.valCrs[ijcrs]*vec.template segment<ndim>(jblk0*ndim); // jblk0!=iblk
    }
  }
}

}

#endif