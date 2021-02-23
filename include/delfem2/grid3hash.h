#ifndef DFM2_GRID3HASH
#define DFM2_GRID3HASH

namespace delfem2 {

class GridCoordinate3 {
public:
  int i, j, k;
public:
  GridCoordinate3(int i_, int j_, int k_) :
      i(i_), j(j_), k(k_) {}

  bool operator==(const GridCoordinate3 &rhs) const {
    return rhs.i == i && rhs.j == j && rhs.k == k;
  }
};

}


namespace std {
template<>
struct hash<delfem2::GridCoordinate3> {
  std::size_t operator()(const delfem2::GridCoordinate3 &k) const {
    /*
    std::uint32_t hi = std::hash<int>()(k.i);
    std::uint32_t hj = std::hash<int>()(k.j) << 1u;
    std::uint32_t hk = std::hash<int>()(k.k) << 1u;
    return ((hi ^ hj) >> 1u) ^ hk;
     */
    return ((hash<int>()(k.i) ^ (hash<int>()(k.j) << 1)) >> 1) ^
           (hash<int>()(k.k) << 1);
  }
};
}

namespace delfem2 {
double WGrid1Cubic(double x) {
  const double x_abs = std::abs(x);
  if (x_abs < 1)
    return std::pow(x_abs, 3) / 2.0 - std::pow(x_abs, 2) + 2.0 / 3.0;
  else if (x_abs < 2)
    return -std::pow(x_abs, 3) / 6.0 + std::pow(x_abs, 2) - 2.0 * x_abs +
           4.0 / 3.0;
  else
    return 0.0;
};

double dWGrid1Cubic(double x) {
  const double x_abs = std::abs(x);
  if (x_abs < 1)
    return x * x_abs * 3.0 / 2.0 - 2.0 * x;
  else if (x_abs < 2)
    return -x * x_abs / 2.0 + 2.0 * x - 2.0 * x / x_abs;
  else
    return 0.0;
}

double WGrid3Cubic(
    double h_inverse,
    double x_offset,
    double y_offset,
    double z_offset)
{
  return WGrid1Cubic(h_inverse * x_offset)
         * WGrid1Cubic(h_inverse * y_offset)
         * WGrid1Cubic(h_inverse * z_offset);
};

double dWGrid3Cubic(
    double h_inverse,
    double d_offset,
    double offset1,
    double offset2)
{
  return dWGrid1Cubic(h_inverse * d_offset)
         * WGrid1Cubic(h_inverse * offset1)
         * WGrid1Cubic(h_inverse * offset2) * h_inverse;
}

}

#endif