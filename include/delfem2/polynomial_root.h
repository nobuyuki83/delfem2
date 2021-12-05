//
// Created by Nobuyuki Umetani on 2021/12/05.
//

#ifndef POLYNOMIAL_ROOT_H_
#define POLYNOMIAL_ROOT_H_

namespace delfem2 {

template<int n>
int StrumNumber(double x, const double s[n][n]) {
  double v[n];
  for (int i = 0; i < n; ++i) {
    v[i] = s[i][n - 1 - i];
    for (int j = 1; j < n - i; ++j) {
      v[i] = v[i] * x + s[i][n - 1 - i - j];
    }
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

int NumberOfRootsInRange_QuadraticPolynomial(
  double x0, double x1,
  double a, double b, double c) {
  assert(x0 <= x1);
  const double sign_a = (a > 0) ? 1. : -1.;
  const double strum[3][3] = {
    {c, b, a},
    {b, 2 * a, 0},
    {sign_a * (b * b - 4 * c * a), 0, 0}
  };
  return StrumNumber<3>(x0, strum) - StrumNumber<3>(x1, strum);
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

}

#endif //POLYNOMIAL_ROOT_H_
