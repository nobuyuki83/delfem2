//
// Created by Nobuyuki Umetani on 2021/12/05.
//

#ifndef POLYNOMIAL_ROOT_H_
#define POLYNOMIAL_ROOT_H_

#include <stack>
#include <tuple>

namespace delfem2 {

/**
 * Evaluate polynomial at x
 * @tparam n degree + 1
 * @param x
 * @param a  coefficient of polynomial f(x)=a[0]+a[1]*x^1+a[2]*a^2 ...
 * @return f(x)
 */
inline double Eval_Polynomial(double x, const double *a, int n) {
  double v = a[n - 1];
  for (int i = 1; i < n; ++i) {
    v = v * x + a[n - 1 - i];
  }
  return v;
}

template<int n>
int StrumNumber(double x, const double s[n][n]) {
  double v[n];
  for (int i = 0; i < n; ++i) {
    v[i] = Eval_Polynomial(x, s[i], n - i);
  }
  // -----------
  int root_number = 0;
  double prev = 0.;
  for (unsigned int i = 0; i < n; ++i) {
    if (v[i] != 0.0) {
      if (prev != 0. && v[i] * prev < 0) { ++root_number; }
      prev = v[i];
    }
  }
  return root_number;
}

/**
 * quadratic polynomial
 * @param strum
 * @param coe
 */
inline void StrumSequenceOfPolynomial_Quadratic(
  double strum[3][3],
  const double coe[3]) {
  const double sign_a = (coe[2] > 0) ? 1. : -1.;
  strum[0][0] = coe[0];
  strum[0][1] = coe[1];
  strum[0][2] = coe[2];
  strum[1][0] = coe[1];
  strum[1][1] = 2 * coe[2];
  strum[2][0] = sign_a * (coe[1] * coe[1] - 4 * coe[0] * coe[2]);
}

/**
 * cubic polynomial
 * @param strum
 * @param coe coe[0] + coe[1]*t + coe[2]*t^2 + coe[3]*t^3
 */
inline void StrumSequenceOfPolynomial_Cubic(
  double strum[4][4],
  const double coe[4]) {
  const double sign_a = (coe[3] > 0) ? 1. : -1.;
  strum[0][0] = coe[0];
  strum[0][1] = coe[1];
  strum[0][2] = coe[2];
  strum[0][3] = coe[3];
  strum[1][0] = coe[1];
  strum[1][1] = 2 * coe[2];
  strum[1][2] = 3 * coe[3];
  strum[2][0] = sign_a * (coe[2] * coe[1] - 9 * coe[3] * coe[0]);
  strum[2][1] = sign_a * (2 * coe[2] * coe[2] - 6 * coe[3] * coe[1]);
  strum[3][0] =
    -coe[3] * (
      -coe[2] * coe[2] * coe[1] * coe[1]
        + 4 * coe[2] * coe[2] * coe[2] * coe[0]
        - 18 * coe[3] * coe[2] * coe[1] * coe[0]
        + coe[3] * (4 * coe[1] * coe[1] * coe[1] + 27 * coe[3] * coe[0] * coe[0])
    );
}

template<int n>
inline void StrumSequenceOfPolynomial(
  double strum[n][n],
  const double coe[n])
{
  static_assert(n>=1, "degree of polynominal should not be negative");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == 0) {
        strum[i][j] = coe[j];
      } else if (i == 1 && j < n - 1) {
        strum[i][j] = (j + 1) * coe[j + 1];
      } else {
        strum[i][j] = 0.;
      }
    }
  }

  for (int i = 0; i < n - 2; i++) {
    int j = i + 1, k1 = n - i - 1, k2 = n - j - 1;
    while (strum[i][k1] == 0) {
      k1--;
      if (k1 < 0) { return; }
    }
    while (strum[j][k2] == 0) {
      k2--;
      if (k2 < 0) { return; }
    }
    double poly[n];
    for (int l = 0; l < n; l++) {
      poly[l] = strum[i][l];
    }
    while (k1 >= k2) {
      while (poly[k1] == 0.) {
        k1--;
        if (k1 < k2) { break; }
      }
      if (k1 >= k2) {
        double quotient = poly[k1] / strum[j][k2];
        poly[k1] = 0.;
        for (int l = 1; l <= k2; l++) {
          poly[k1 - l] -= quotient * strum[j][k2 - l];
        }
      }
    }
    for (int l = 0; l < n; l++) {
      strum[j + 1][l] = -poly[l];
    }
  }
}

template<int n>
inline std::vector<std::pair<double, double> > RootInterval_StrumSequence(
  double x0, double x1,
  const double strum[n][n]) {

  std::stack<std::tuple<double, double, int, int> > next;
  next.emplace(
    x0, x1,
    delfem2::StrumNumber<n>(x0, strum),
    delfem2::StrumNumber<n>(x1, strum));

  std::vector<std::pair<double, double>> res;
  while (!next.empty()) {
    const auto intvl = next.top();
    next.pop();
    const int n_left = std::get<2>(intvl);
    const int n_right = std::get<3>(intvl);
    const double x_left = std::get<0>(intvl);
    const double x_right = std::get<1>(intvl);
    if (n_left - n_right > 1) {  // do the bi-section
      const double x_mid = (x_left + x_right) / 2;
      const int n_mid = delfem2::StrumNumber<n>(x_mid, strum);
      if (n_left > n_mid) {  // left half interval
        next.emplace(x_left, x_mid, n_left, n_mid);
      }
      if (n_mid > n_right) {  // right half interval
        next.emplace(x_mid, x_right, n_mid, n_right);
      }
    } else if (n_left - n_right == 1) {
      res.emplace_back(x_left, x_right);
    }
  }
  return res;
}

template<typename REAL>
std::vector<REAL> RootsInRange_QuadraticPolynomial(
  REAL x0, REAL x1,
  REAL a, REAL b, REAL c) {
  assert(x0 <= x1);
  REAL delta = b * b - 4 * a * c;
  if (delta < 0) {
    return {};
  }
  if (a == 0) {
    if (b == 0) {
      assert(0);
      return {};
    }
    const REAL cand = c / b;
    if (cand > x0 || cand < x1) { return {cand}; }
    return {};
  }
  delta = std::sqrt(delta);
  const REAL cands[2] = {
    (-b - delta) / (2 * a),
    (-b + delta) / (2 * a)};
  std::vector<REAL> res;
  res.reserve(2);
  for (auto cand: cands) {
    if (cand < x0 || cand > x1) { continue; }
    res.push_back(cand);
  }
  return res;
}

/**
 *
 * @param range_left  lower bound of range
 * @param range_right  upper bound of range
 * @param a coefficient of polynomial
 * @param ndegree_plus1 degree of polynomial plus 1
 * @param num_iteration  number of iteration
 * @return
 */
inline double RootInInterval_Bisection(
  double range_left,
  double range_right,
  const double *a,
  int ndegree_plus1,
  int num_iteration) {
  double x0 = range_left;
  double x1 = range_right;
  double fx0 = Eval_Polynomial(x0, a, ndegree_plus1);
  double fx1 = Eval_Polynomial(x1, a, ndegree_plus1);
  assert(fx0 * fx1 < 0);
  for (int i = 0; i < num_iteration; i++) {
    const double x2 = (x0 + x1) / 2; // bisection method
    const double fx2 = Eval_Polynomial(x2, a, ndegree_plus1);
    if (fx1 * fx2 < 0) {
      x0 = x1;
      x1 = x2;
      fx0 = fx1;
      fx1 = fx2;
    } else {
      x1 = x2;
      fx1 = fx2;
    }
  }
  return (x0 * fx1 - x1 * fx0) / (fx1 - fx0); // secant method
}

 /**
  * Find the root of polynomial in the range [0,1]
  * @tparam[in] ndegplus1 (degree of polynomial) + 1, e.g.(ndegplus1==3)->quadratic
  * @param[in] coe f(t) = coe[0] + coe[1]*t + coe[2]*t^2 + ...
  * @param[in] num_bisection
  * @return roots of polynomial
  */
template<int ndegplus1>
std::vector<double> RootsOfPolynomial(
  const double coe[ndegplus1],
  int num_bisection){
  double strum[ndegplus1][ndegplus1];
  StrumSequenceOfPolynomial(strum, coe);
  std::vector<std::pair<double, double>> intvls = RootInterval_StrumSequence<ndegplus1>(0, 1, strum);
  std::vector<double> res;
  for (auto intvl: intvls) {
    double y0 = RootInInterval_Bisection(
      intvl.first, intvl.second,
      coe, ndegplus1, num_bisection);
    res.push_back(y0);
  }
  return res;
}


}

#endif //POLYNOMIAL_ROOT_H_
