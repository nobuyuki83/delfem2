#ifndef DFM2_CUDA_RANDOM_H
#define DFM2_CUDA_RANDOM_H


namespace delfem2 {
namespace cuda {


__device__
uint4 device_Xorshift128(uint4 seed) {
  const uint t = (seed.x ^ (seed.x << 11));
  uint4 res;
  res.x = seed.y;
  res.y = seed.z;
  res.z = seed.w;
  res.w = (seed.w ^ (seed.w >> 19)) ^ (t ^ (t >> 8));
  return res;
}

__device__
float device_RandomRange(uint iker, float min, float max) {
  // 4294967295 is the maximum of 4bytes
  return min + (max - min) * (float(iker) / float(4294967295));
}


}
}

#endif