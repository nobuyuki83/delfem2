//
// Created by Nobuyuki Umetani on 2021/12/05.
//

#ifndef POLYNOMIAL_ROOT_H_
#define POLYNOMIAL_ROOT_H_

#include <stack>
#include <tuple>

namespace delfem2 {

/**
 *
 * @tparam n degree + 1
 * @param x
 * @param a  f(x)=a[0]+a[1]*x^1+a[2]*a^2 ...
 * @return f(x)
 */
template<int n>
double Eval_Polynomial(double x, const double a[n]){
  double v = a[n-1];
  for(int i=1;i<n;++i){
    v += (v + a[n-1-i])*x;
  }
  return v;
}

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

void StrumSequence_QuadraticPolynomial(
  double strum[3][3],
  double a, double b, double c) {
  const double sign_a = (a > 0) ? 1. : -1.;
  strum[0][0] = c;
  strum[0][1] = b;
  strum[0][2] = a;
  strum[1][0] = b;
  strum[1][1] = 2 * a;
  strum[2][0] = sign_a * (b * b - 4 * c * a);
}

void StrumSequence_CubicPolynomial(
  double strum[4][4],
  double a, double b, double c, double d) {
  const double sign_a = (a > 0) ? 1. : -1.;
  strum[0][0] = d;
  strum[0][1] = c;
  strum[0][2] = b;
  strum[0][3] = a;
  strum[1][0] = c;
  strum[1][1] = 2*b;
  strum[1][2] = 3*a;
  strum[2][0] = sign_a * (b * c - 9 * a * d);
  strum[2][1] = sign_a * (2 * b * b - 6 * a * c);
  strum[3][0] = - a * (-b*b*c*c + 4*b*b*b*d - 18*a*b*c*d + a * (4*c*c*c + 27*a*d*d));
}

template <int n>
std::vector< std::pair<double, double> > RootInterval_StrumSequence(
  double x0, double x1,
  double strum[n][n]) {

  std::stack< std::tuple<double, double, int, int> > next;
  next.emplace(
    x0,x1,
    delfem2::StrumNumber<n>(x0,strum),
    delfem2::StrumNumber<n>(x1,strum) );

  std::vector<std::pair<double, double>> res;
  while (!next.empty()) {
    const auto intvl = next.top();
    next.pop();
    int nl = std::get<2>(intvl);
    int nr = std::get<3>(intvl);
    if ( nl-nr > 1) {  // do the bi-section
      double xl = std::get<0>(intvl);
      double xr = std::get<1>(intvl);
      int mid_strum = delfem2::StrumNumber<n>((xl + xr) / 2, strum);
      if (nl > mid_strum) {  // left half interval
        next.emplace(xl, (xl + xr) / 2, nl, mid_strum);
      }
      if (mid_strum > nr) {  // right half interval
        next.emplace((xl + xr) / 2, xr, mid_strum, nr);
      }
    } else if (nl-nr == 1) {
      res.emplace_back(std::get<0>(intvl), std::get<1>(intvl));
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
    if( cand > x0 || cand < x1 ){ return {cand}; }
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
