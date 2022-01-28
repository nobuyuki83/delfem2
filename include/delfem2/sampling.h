//
// Created by Nobuyuki Umetani on 2022/01/28.
//

#ifndef DFM2_SAMPLING_H_
#define DFM2_SAMPLING_H_

#include <array>
#include <random>

namespace delfem2 {

// based on https://github.com/postgres/postgres/blob/master/src/port/erand48.c
template<typename T>
T my_erand48(std::array<unsigned short, 3> &xseed) {
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
std::array<T, 3> RandomVec3(
    std::uniform_real_distribution<T> &dist,
    std::mt19937 &reng) {
  return {dist(reng), dist(reng), dist(reng)};
}

template<typename T>
T SampleTent(
    std::array<unsigned short, 3> &Xi) {
  const T r1 = 2 * my_erand48<T>(Xi);
  return r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);  // tent filter (-1 .. +1 )
}

// Sample an hemispherical direction with uniform distribution.
template<typename T>
inline std::array<T, 3> SampleHemisphereZup(
    std::array<unsigned short, 3> &Xi) {
  auto z = my_erand48<T>(Xi);
  auto r = std::sqrt(std::clamp(1 - z * z, 0, 1));
  auto phi = 2 * M_PI * my_erand48<T>(Xi);
  return {r * cos(phi), r * sin(phi), z};
}

// Sample an hemispherical direction with cosine distribution.
template<typename T>
inline std::array<T, 3> SampleHemisphereZupCos(
    std::array<unsigned short, 3> &Xi) {
  auto z = std::sqrt(my_erand48<T>(Xi));
  auto r = std::sqrt(1 - z * z);
  auto phi = 2 * M_PI * my_erand48<T>(Xi);
  return {r * cos(phi), r * sin(phi), z};
}

}

#endif //DFM2_SAMPLING_H_
