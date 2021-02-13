#ifndef DFM2_EIGEN_GRID2_H
#define DFM2_EIGEN_GRID2_H

namespace delfem2 {
namespace eigen {

namespace grid2 {
double Length2D(
    const double *p,
    const double *q) {
  return sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]));
}
}

void Interpolate_Grid_ThinPlateSpline(
    std::vector<float> &aGridZ,
    unsigned int nw,
    unsigned int nh,
    const std::vector<double> &aPntXY,
    const std::vector<float> &aPntZ) {
  const size_t np = aPntXY.size() / 2;
  Eigen::MatrixXd A(np, np);
  A.setZero();
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int jp = 0; jp < np; ++jp) {
      double r = grid2::Length2D(aPntXY.data() + ip * 2, aPntXY.data() + jp * 2);
      A(ip, jp) = r * r * log(r + 1.0e-8);
    }
  }
  Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
  Eigen::VectorXd b(np);
  for (unsigned int ip = 0; ip < np; ++ip) { b(ip) = aPntZ[ip]; }
  Eigen::VectorXd w = lu.solve(b);

  aGridZ.resize(nw * nh);
  for (unsigned int iw = 0; iw < nw; ++iw) {
    for (unsigned int ih = 0; ih < nh; ++ih) {
      const double p0[2] = {double(iw), double(ih)};
      double val = 0.;
      for (unsigned int ip = 0; ip < np; ++ip) {
        double r = grid2::Length2D(p0, aPntXY.data() + ip * 2);
        val += r * r * log(r + 1.0e-8) * w[ip];
      }
      aGridZ[ih * nw + iw] = val;
//      std::cout << ih << " " << iw << " " << val << std::endl;
    }
  }
}

}
}

#endif