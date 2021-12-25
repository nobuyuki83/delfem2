//
// Created by Nobuyuki Umetani on 2021-10-17.
//

#ifndef IMPLICIT_RBF_APPROXIMATION_H_
#define IMPLICIT_RBF_APPROXIMATION_H_

#include <functional>
#include <tuple>
#include <vector>
#include <Eigen/LU>

class ImplicitRbfApproximation {
 public:
  ImplicitRbfApproximation(
      std::function<double(double)> rbf,
      bool is_linear)
      : rbf(std::move(rbf)), is_linear(is_linear) {}

  void Clear(){
    sample_xy.clear();
    weights.clear();
  }

  [[nodiscard]] bool IsInitailized2() const {
    const auto num_sample = static_cast<unsigned int>(sample_xy.size() / 2);
    if( num_sample == 0 ){ return false; }
    const unsigned int ndof = is_linear ? num_sample + 3 : num_sample;
    if(weights.size() != ndof){ return false; }
    return true;
  }

  template<class V2>
  void SetPolyline2(const std::vector<V2> &stroke, double eps) {
    GenerateSamplePositionsForPolyline(stroke, eps);
    const auto num_sample = static_cast<unsigned int>(sample_xy.size() / 2);
    const unsigned int ndof = is_linear ? num_sample + 3 : num_sample;
    Eigen::MatrixXd A(ndof, ndof);
    Eigen::VectorXd y(ndof);
    A.setZero();
    y.setZero();
    for (unsigned int i = 0; i < num_sample; i++) {
      for (unsigned int j = 0; j < num_sample; j++) {
        const V2 pi(sample_xy[i * 2 + 0], sample_xy[i * 2 + 1]);
        const V2 pj(sample_xy[j * 2 + 0], sample_xy[j * 2 + 1]);
        double r = (pi - pj).norm();
        A(i, j) = rbf(r);
      }
      y(i) = (i % 2 == 0) ? eps : -eps;
    }
    if (is_linear) {
      for (unsigned int ismpl = 0; ismpl < num_sample; ismpl++) {
        A(num_sample + 0, ismpl) = 1;
        A(num_sample + 1, ismpl) = sample_xy[ismpl * 2 + 0];
        A(num_sample + 2, ismpl) = sample_xy[ismpl * 2 + 1];
        A(ismpl, num_sample + 0) = 1;
        A(ismpl, num_sample + 1) = sample_xy[ismpl * 2 + 0];
        A(ismpl, num_sample + 2) = sample_xy[ismpl * 2 + 1];
      }
    }
    Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
    Eigen::VectorXd x = solver.solve(y);
    weights.resize(ndof);
    for (unsigned int idof = 0; idof < ndof; idof++) {
      weights[idof] = x(idof);
    }
  }

  [[nodiscard]] double Evaluate2(double x, double y) const {
    const size_t num_samples = sample_xy.size() / 2;
    const unsigned int ndof = is_linear ? num_samples + 3 : num_samples;
    assert(weights.size() == ndof);
    double t = 0;
    for (unsigned int ismpl = 0; ismpl < num_samples; ismpl++) {
      double dx0 = sample_xy[ismpl * 2 + 0] - x;
      double dy0 = sample_xy[ismpl * 2 + 1] - y;
      double r = std::sqrt(dx0 * dx0 + dy0 * dy0);
      t += weights[ismpl] * rbf(r);
    }
    if (is_linear) {
      t += weights[num_samples] + weights[num_samples + 1] * x + weights[num_samples + 2] * y;
    }
    return t;
  }

  [[nodiscard]] std::tuple<double,double,double> Evaluate2Grad(
      double x, double y,
      const std::function<double(double)> &diff_rbf) const {
    const size_t num_samples = sample_xy.size() / 2;
    const unsigned int ndof = is_linear ? num_samples + 3 : num_samples;
    assert(weights.size() == ndof);
    double t = 0., dtdx = 0., dtdy = 0.;
    for (unsigned int ismpl = 0; ismpl < num_samples; ismpl++) {
      const double dx0 = x - sample_xy[ismpl * 2 + 0];
      const double dy0 = y - sample_xy[ismpl * 2 + 1];
      const double r = std::sqrt(dx0 * dx0 + dy0 * dy0);
      t += weights[ismpl] * rbf(r);
      dtdx += weights[ismpl] * diff_rbf(r) * dx0 / r;
      dtdy += weights[ismpl] * diff_rbf(r) * dy0 / r;
    }
    if (is_linear) {
      t += weights[num_samples] + weights[num_samples + 1] * x + weights[num_samples + 2] * y;
      dtdx += weights[num_samples + 1];
      dtdy += weights[num_samples + 2];
    }
    return {t,dtdx,dtdy};
  }

 private:
  template<class V2>
  void GenerateSamplePositionsForPolyline(const std::vector<V2> &stroke, double eps) {
    const size_t n = stroke.size();
    const size_t num_sample = n * 2;
    std::vector<V2> aNorm;
    aNorm.resize(n, V2(0, 0));
    for (unsigned int is = 0; is < stroke.size() - 1; is++) {
      V2 p0 = stroke[is + 0];
      V2 p1 = stroke[is + 1];
      V2 e = (p1 - p0).normalized();
      V2 v(e.y, -e.x);
      aNorm[is + 0] += v;
      aNorm[is + 1] += v;
    }
    for (unsigned int i = 0; i < n; i++) { aNorm[i].normalize(); }
    sample_xy.resize(num_sample * 2);
    for (unsigned int ismpl = 0; ismpl < num_sample; ismpl++) {
      unsigned int i0 = ismpl / 2;
      V2 pi = stroke[i0];
      pi += (ismpl % 2 == 0) ? aNorm[i0] * eps : -1.0 * aNorm[i0] * eps;
      sample_xy[ismpl * 2 + 0] = pi[0];
      sample_xy[ismpl * 2 + 1] = pi[1];
    }
  }
 private:
  const std::function<double(double)> rbf;
  const bool is_linear;
  //
  std::vector<double> sample_xy;
  std::vector<double> weights;
};


#endif //IMPLICIT_RBF_APPROXIMATION_H_
