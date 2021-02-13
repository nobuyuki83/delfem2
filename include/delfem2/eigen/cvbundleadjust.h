#ifndef DFM2_EIGEN_CVBUNDLEADJUST
#define DFM2_EIGEN_CVBUNDLEADJUST

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/cvbundleadjust.h"
#include <vector>
#include <Eigen/Core>
#include <Eigen/LU>

namespace delfem2 {
namespace eigen {

double BundleAdjust_ImgPair_NewtonIter(
    std::vector<float> &aZ0,
    double *R01,
    double *t01,
    //
    const std::vector<double> &aPnt0,
    const double *K0inv,
    const double *K1,
    const std::vector<double> &aPnt1,
    double eps) {
  namespace dfm2 = delfem2;
  const size_t np = aPnt0.size() / 2;
  assert(aPnt1.size() == np * 2);
  //
  std::vector<double> F, dF;
  double E = ReprojectionEnergy(
      F, dF,
      aPnt0, K0inv, aZ0, R01, t01, K1, aPnt1);
  const size_t M = np + 6;
  assert(F.size() == np * 2);
  assert(dF.size() == np * 2 * 7);
  Eigen::MatrixXd A(M, M);
  A.setZero();
  for (unsigned int ip = 0; ip < np; ++ip) {
    A(ip, ip) = dF[ip * 14 + 0] * dF[ip * 14 + 0] + dF[ip * 14 + 7] * dF[ip * 14 + 7];
    for (int j = 0; j < 6; ++j) {
      A(ip, np + j) = dF[ip * 14 + 0] * dF[ip * 14 + 1 + j] + dF[ip * 14 + 7] * dF[ip * 14 + 8 + j];
    }
  }
  for (int i = 0; i < 6; ++i) {
    for (unsigned int jp = 0; jp < np; ++jp) {
      A(np + i, jp) = dF[jp * 14 + 1 + i] * dF[jp * 14 + 0] + dF[jp * 14 + 8 + i] * dF[jp * 14 + 7];
    }
    for (int j = 0; j < 6; ++j) {
      for (unsigned int kp = 0; kp < np; ++kp) {
        A(np + i, np + j) += dF[kp * 14 + 1 + i] * dF[kp * 14 + 1 + j] + dF[kp * 14 + 8 + i] * dF[kp * 14 + 8 + j];
      }
    }
  }
  for (unsigned int i = 0; i < M; ++i) { A(i, i) += eps; }
  Eigen::VectorXd b = Eigen::VectorXd::Zero(M);
  for (unsigned int ip = 0; ip < np; ++ip) {
    b(ip) = F[ip * 2 + 0] * dF[ip * 14 + 0] + F[ip * 2 + 1] * dF[ip * 14 + 7];
    for (int i = 0; i < 6; ++i) {
      b(np + i) += F[ip * 2 + 0] * dF[ip * 14 + 1 + i] + F[ip * 2 + 1] * dF[ip * 14 + 8 + i];
    }
  }
  b *= -1.0;
  Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
  Eigen::VectorXd x = lu.solve(b);
  //std::cout << E << std::endl;
  for (unsigned int ip = 0; ip < np; ++ip) {
    aZ0[ip] += x(ip);
  }
  t01[0] += x(np + 0);
  t01[1] += x(np + 1);
  t01[2] += x(np + 2);
  double dw[3] = {x(np + 3), x(np + 4), x(np + 5)};
  double dR[9];
  dfm2::Mat3_Rotation_Cartesian(dR, dw);
  double R0[9];
  dfm2::Copy_Mat3(R0, R01);
  dfm2::MatMat3(R01, dR, R0);
  return E;
}

}
}

#endif