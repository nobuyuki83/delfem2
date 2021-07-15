
#include <vector>
#include <Eigen/Core>

namespace delfem2 {

// {y_vec} = \beta * {y_vec} + \alpha * [mat] * {p_vec}
void AddMatVec(
    Eigen::VectorXd &y_vec,
    double beta,
    double alpha,
    const Eigen::MatrixXd &mat,
    const Eigen::VectorXd &x_vec) {
  y_vec = alpha * mat * x_vec + beta * y_vec;
}


template <int nrow, int ncol, int ndimrow, int ndimcol, typename T>
bool Merge(
    Eigen::MatrixXd& A,
    const unsigned int* aIpRow,
    const unsigned int* aIpCol,
    const T emat[nrow][ncol][ndimrow][ndimcol],
    std::vector<unsigned int>& tmp_buffer)
{
  for(unsigned int in=0;in<nrow;++in){
    for(unsigned int jn=0;jn<ncol;++jn) {
      const unsigned int ip = aIpRow[in];
      const unsigned int jp = aIpCol[jn];
      for(unsigned int idim=0;idim<ndimrow;++idim){
        for(unsigned int jdim=0;jdim<ndimcol;++jdim) {
          A(ip*ndimrow+idim,jp*ndimcol+jdim) += emat[in][jn][idim][jdim];
        }
      }
    }
  }
  return true;
}

// above: Eigen::MatrixX
// ------------------------------
// below: Eigen::VectorX


template <typename REAL>
void AddScaledVec(
    Eigen::Matrix<REAL,-1,1> &y,
    REAL alpha,
    const Eigen::Matrix<REAL,-1,1> &x) {
  assert(y.rows() == x.rows());
  y += alpha * x;
}

template<typename REAL, int nDim>
void AddScaledVec(
    Eigen::Matrix<REAL,-1,nDim,Eigen::RowMajor,-1,nDim> &y,
    REAL alpha,
    const Eigen::Matrix<REAL,-1,nDim,Eigen::RowMajor,-1,nDim> &x) {
  assert(y.rows() == x.rows());
  y += alpha * x;
}

template<typename REAL>
void ScaleAndAddVec(
    Eigen::Matrix<REAL,-1,1> &y,
    REAL beta,
    const Eigen::Matrix<REAL,-1,1> &x) {
  assert(y.cols() == x.cols());
  assert(y.rows() == x.rows());
  y = beta * y + x;
}


template<typename REAL, int nDim>
void ScaleAndAddVec(
    Eigen::Matrix<REAL,-1,nDim,Eigen::RowMajor,-1,nDim> &y,
    REAL beta,
    const Eigen::Matrix<REAL,-1,nDim,Eigen::RowMajor,-1,nDim> &x) {
  assert(y.rows() == x.rows());
  y = beta * y + x;
}

// --------------------------------------

template<typename REAL>
REAL Dot(
    const Eigen::Matrix<REAL,-1,1> &y,
    const Eigen::Matrix<REAL,-1,1> &x) {
  assert(y.rows() == x.rows());
  return y.dot(x);
}


template<typename REAL, int nDim>
REAL Dot(
    const Eigen::Matrix<REAL,-1,nDim,Eigen::RowMajor,-1,nDim> &y,
    const Eigen::Matrix<REAL,-1,nDim,Eigen::RowMajor,-1,nDim> &x){
  assert(y.rows() == x.rows());
  return y.cwiseProduct(x).sum();
}

// --------------------------------------------

template<class VEC>
void setZero_Flag(
    VEC &vec_b,
    const std::vector<int> &aBCFlag,
    int iflag_nonzero) {
  const std::size_t ndof = vec_b.size();
  for (unsigned int i = 0; i < ndof; ++i) {
    if (aBCFlag[i] == iflag_nonzero) continue;
    vec_b(i) = 0;
  }
}

template<typename REAL, int nDim>
void setZero_Flag(
    Eigen::Matrix<REAL,-1,nDim,Eigen::RowMajor,-1,nDim> &vec_b,
    unsigned int np,
    const std::vector<int> &aBCFlag,
    int iflag_nonzero)
{
  const std::size_t ndim = aBCFlag.size()/np;
  assert( vec_b.size() % np == 0 );
  assert( ndim <= nDim );
  for (unsigned int ip = 0; ip < np; ++ip) {
    for(unsigned int idim=0;idim<ndim;++idim) {
      if (aBCFlag[ip * ndim + idim] == iflag_nonzero) continue;
      vec_b(ip,idim) = 0;
    }
  }
}

// -------------------------------------------

template<typename REAL>
void XPlusAY(
    std::vector<REAL> &X,
    const std::vector<int> &aBCFlag,
    REAL alpha,
    const Eigen::Matrix<REAL,-1,1> &Y)
{
  const unsigned int nDoF = aBCFlag.size();
  assert(nDoF == Y.rows());
  for (unsigned int i = 0; i < nDoF; ++i) {
    if (aBCFlag[i] != 0) continue;
    X[i] += alpha * Y[i];
  }
}

template<typename REAL, int nDim>
void XPlusAY(
    std::vector<REAL> &X,
    const unsigned int np,
    const std::vector<int> &aBCFlag,
    REAL alpha,
    Eigen::Matrix<REAL,-1,nDim,Eigen::RowMajor,-1,nDim> &Y)
{
  const std::size_t ndim = aBCFlag.size()/np;
  assert( aBCFlag.size() % np == 0 );
  assert(X.size() == aBCFlag.size() );
  assert( ndim <= nDim );
  for (unsigned int ip = 0; ip < np; ++ip) {
    for(unsigned int idim=0;idim<ndim;++idim){
      if (aBCFlag[ip*ndim+idim] != 0) continue;
      X[ip*ndim+idim] += alpha * Y(ip,idim);
    }
  }
}

}

//#include "delfem2/lsitrsol.h"
