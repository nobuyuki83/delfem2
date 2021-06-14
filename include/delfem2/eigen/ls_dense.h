
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

template<typename REAL>
void AddScaledVec(
    Eigen::VectorX<REAL> &y,
    REAL alpha,
    const Eigen::VectorX<REAL> &x) {
  assert(y.cols() == x.cols());
  y += alpha * x;
}

template<typename REAL>
void ScaleAndAddVec(
    Eigen::VectorX<REAL> &y,
    REAL beta,
    const Eigen::VectorX<REAL> &x) {
  assert(y.cols() == x.cols());
  y = beta * y + x;
}


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

template<typename REAL>
void XPlusAY(
    std::vector<REAL> &X,
    const unsigned int nDoF,
    const std::vector<int> &aBCFlag,
    REAL alpha,
    const Eigen::VectorX<REAL> &Y) {
  for (unsigned int i = 0; i < nDoF; ++i) {
    if (aBCFlag[i] != 0) continue;
    X[i] += alpha * Y[i];
  }
}

}

//#include "delfem2/lsitrsol.h"
