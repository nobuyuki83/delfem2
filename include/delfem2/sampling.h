//
// Created by Nobuyuki Umetani on 2022/01/28.
//

#ifndef DFM2_SAMPLING_H_
#define DFM2_SAMPLING_H_

#include <array>
#include <random>
#include <algorithm>  // clamp

#ifndef M_PI
#  define M_PI 3.1415926535897932384626433832
#endif

namespace delfem2 {

// based on https://github.com/postgres/postgres/blob/master/src/port/erand48.c
template<typename T>
T MyERand48(std::array<unsigned short, 3> &xseed) {
  constexpr unsigned short my_rand48_mult[3] = {
      0xe66d,
      0xdeec,
      0x0005
  };
  constexpr unsigned short my_rand48_add = 0x000b;

  unsigned long accu;
  unsigned short temp[2];

  accu = (unsigned long) my_rand48_mult[0] * (unsigned long) xseed[0] +
      (unsigned long) my_rand48_add;
  temp[0] = (unsigned short) accu;    /* lower 16 bits */
  accu >>= sizeof(unsigned short) * 8;
  accu += (unsigned long) my_rand48_mult[0] * (unsigned long) xseed[1] +
      (unsigned long) my_rand48_mult[1] * (unsigned long) xseed[0];
  temp[1] = (unsigned short) accu;    /* middle 16 bits */
  accu >>= sizeof(unsigned short) * 8;
  accu += my_rand48_mult[0] * xseed[2] + my_rand48_mult[1] * xseed[1] + my_rand48_mult[2] * xseed[0];
  xseed[0] = temp[0];
  xseed[1] = temp[1];
  xseed[2] = (unsigned short) accu;
  // --------
  return
      std::ldexp(static_cast<T>(xseed[0]), -48) +
          std::ldexp(static_cast<T>(xseed[1]), -32) +
          std::ldexp(static_cast<T>(xseed[2]), -16);
}

template<typename T>
std::array<T, 2> RandomVec2(
    std::array<unsigned short, 3> &xi) {
  return { MyERand48<T>(xi), MyERand48<T>(xi) };
}

template<int n, typename T>
std::array<T, n> RandomVec(
    std::uniform_real_distribution<T> &dist,
    std::mt19937 &reng,
    T mag = 1.0,
    T offset = 0.0) {
  std::array<T, n> r;
  for(int i=0;i<n;++i){
    r[i] =mag * dist(reng) + offset;
  }
  return r;
}

template<typename T>
std::array<T, 9> RandomMat3(
    std::uniform_real_distribution<T> &dist,
    std::mt19937 &reng,
    T mag = 1.0,
    T offset = 0.0) {
  return {
      mag * dist(reng) + offset,
      mag * dist(reng) + offset,
      mag * dist(reng) + offset,
      mag * dist(reng) + offset,
      mag * dist(reng) + offset,
      mag * dist(reng) + offset,
      mag * dist(reng) + offset,
      mag * dist(reng) + offset,
      mag * dist(reng) + offset };
}

template<typename T>
T SampleTent(
    std::array<unsigned short, 3> &Xi) {
  const T r1 = 2 * MyERand48<T>(Xi);
  return r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);  // tent filter (-1 .. +1 )
}

// Sample an hemispherical direction with uniform distribution.
template<typename T>
inline std::array<T, 3> SampleHemisphereZupUniform(
    const std::array<T, 2> &v2) {
  auto z = v2[0];
  auto r = std::sqrt(std::clamp(1 - z * z, 0, 1));
  auto phi = 2 * M_PI * v2[1];
  return {r * cos(phi), r * sin(phi), z};
}

// Sample an hemispherical direction with cosine distribution.
template<typename T>
inline std::array<T, 3> SampleHemisphereZupCos(
    const std::array<T, 2> &v2) {
  auto z = std::sqrt(v2[0]);
  auto r = std::sqrt(1 - z * z);
  auto phi = 2 * M_PI * v2[1];
  return {r * cos(phi), r * sin(phi), z};
}

// Sample an hemispherical direction with cosine distribution.
template<typename VEC, typename SCALAR = typename VEC::Scalar>
inline VEC SampleHemisphereNormalCos(
    const VEC &n,
    const std::array<SCALAR, 2> &v2) {
  const std::array<SCALAR, 3> h0 = SampleHemisphereZupCos(v2);
  const VEC u = ((fabs(n[0]) > .1 ? VEC(0, 1, 0) : VEC(1,0,0)).cross(n)).normalized();  // orthogonal to w
  const VEC v = n.cross(u);  // orthogonal to w and u
  return (h0[0] * u + h0[1] * v + h0[2] * n).normalized();
}

}

#endif //DFM2_SAMPLING_H_
