//
// Created by Nobuyuki Umetani on 2022/02/10.
//

#ifndef DFM2_FDM_ARRAY2_H_
#define DFM2_FDM_ARRAY2_H_

#include <vector>
#include <array>

template<typename T>
class FdmArray2 {
 public:
  FdmArray2(int ni_, int nj_) : ni(ni_), nj(nj_) { v.resize(ni * nj); }
  FdmArray2(int ni_, int nj_, T v_) : ni(ni_), nj(nj_) { v.resize(ni * nj, v_); }

  const T &operator()(int i, int j) const {
    assert(i >= 0 && i < ni && j >= 0 && j < nj);
    return v[i + ni * j];
  }

  T &operator()(int i, int j) {
    assert(i >= 0 && i < ni && j >= 0 && j < nj);
    return v[i + ni * j];
  }

  // Clamped Fetch
  T ClampedFetch(
      int i,
      int j) const {
    i = myclamp(i, 0, ni-1);
    j = myclamp(j, 0, nj-1);
    assert(i >= 0 && i < ni && j >= 0 && j < nj);
    return v[ni * j + i];
  }

  T& ClampedFetch(
      int i,
      int j) {
    i = myclamp(i, 0, ni-1);
    j = myclamp(j, 0, nj-1);
    assert(i >= 0 && i < ni && j >= 0 && j < nj);
    return v[ni * j + i];
  }
 private:
  template <typename S>
  static S mymax(S i, S j) { return i > j ? i : j; };

  template <typename S>
  static S mymin(S i, S j) { return i < j ? i : j; };

  template <typename S>
  static S myclamp(S iv, S imin, S imax){ return mymin(mymax(imin, iv), imax); }
 public:
  int ni, nj;
  std::vector<T> v;
};

template<int ndim>
double LinearInterpolationOnGrid2(
    const FdmArray2<std::array<double, ndim>> &data,
    unsigned int idim,
    double x,
    double y) {
  unsigned int ni = data.ni;
  unsigned int nj = data.nj;
  auto mymin = [](double a, double b) { return a < b ? a : b; };
  auto mymax = [](double a, double b) { return a > b ? a : b; };
  double x1 = mymax(0.0, mymin((double) ni-1.0-1.0e-10, x));
  double y1 = mymax(0.0, mymin((double) nj-1.0-1.0e-10, y));
  auto i = static_cast<unsigned int>(x1);
  auto j = static_cast<unsigned int>(y1);
  assert(i >= 0 && i < (unsigned int)data.ni - 1);
  assert(j >= 0 && j < (unsigned int)data.nj - 1);
  double v00 = data(i + 0, j + 0)[idim];
  double v10 = data(i + 1, j + 0)[idim];
  double v01 = data(i + 0, j + 1)[idim];
  double v11 = data(i + 1, j + 1)[idim];
  double rx = x - i;
  double ry = y - j;
  return (1 - rx) * (1 - ry) * v00 + rx * (1 - ry) * v10 + (1 - rx) * ry * v01 + rx * ry * v11;
}

#endif //DFM2_FDM_ARRAY2_H_
