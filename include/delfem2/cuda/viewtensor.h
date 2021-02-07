#ifndef DFM2_CUDA_VIEWTENSOR_H
#define DFM2_CUDA_VIEWTENSOR_H

namespace delfem2 {
namespace cuda {

template<typename T>
class ViewTensor3 {
public:
  __device__
  ViewTensor3(T *ptr_, unsigned int n0_, unsigned int n1_, unsigned int n2_)
      : ptr(ptr_), s0(n1_ * n2_), s1(n2_) {}

  __device__
  const T &operator()(unsigned int i0, unsigned int i1, unsigned int i2) const {
    return ptr[i0 * s0 + i1 * s1 + i2];
  }

  __device__
  T &operator()(unsigned int i0, unsigned int i1, unsigned int i2) {
    return ptr[i0 * s0 + i1 * s1 + i2];
  }

public:
  const unsigned int s0 = 0, s1 = 0;
  T *const ptr;
};

template<typename T>
class ViewTensor3Const {
public:
  __device__
  ViewTensor3Const(const T *ptr_, unsigned int n0_, unsigned int n1_, unsigned int n2_)
      : ptr(ptr_), s0(n1_ * n2_), s1(n2_) {}

  __device__
  const T &operator()(unsigned int i0, unsigned int i1, unsigned int i2) const {
    return ptr[i0 * s0 + i1 * s1 + i2];
  }

public:
  const unsigned int s0 = 0, s1 = 0;
  const T *const ptr;
};

template<typename T>
class ViewTensor4 {
public:
  __device__
  ViewTensor4(T *ptr_, unsigned int n0_, unsigned int n1_, unsigned int n2_, unsigned int n3_)
      : ptr(ptr_), s0(n1_ * n2_ * n3_), s1(n2_ * n3_), s2(n3_) {}

  __device__
  T *data(unsigned int i0, unsigned int i1, unsigned int i2) {
    return ptr + i0 * s0 + i1 * s1 + i2 * s2;
  }

  __device__
  const T *data(unsigned int i0, unsigned int i1, unsigned int i2) const {
    return ptr + i0 * s0 + i1 * s1 + i2 * s2;
  }

  __device__
  T &operator()(unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3) {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3];
  }

  const T &operator()(unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3) const {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3];
  }

public:
  const unsigned int s0 = 0, s1 = 0, s2 = 0;
  T *const ptr;
};

template<typename T>
class ViewTensor4Const {
public:
  __device__
  ViewTensor4Const(
      const T *ptr_,
      unsigned int n0_,
      unsigned int n1_,
      unsigned int n2_,
      unsigned int n3_)
      : ptr(ptr_), s0(n1_ * n2_ * n3_), s1(n2_ * n3_), s2(n3_) {}

  __device__
  const T *data(
      unsigned int i0,
      unsigned int i1,
      unsigned int i2) const {
    return ptr + i0 * s0 + i1 * s1 + i2 * s2;
  }

  __device__
  const T &operator()(
      unsigned int i0,
      unsigned int i1,
      unsigned int i2,
      unsigned int i3) const {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3];
  }

public:
  const unsigned int s0 = 0, s1 = 0, s2 = 0;
  const T *const ptr;
};

template<typename T>
class ViewTensor5 {
public:
  __device__
  ViewTensor5(
      T *ptr_,
      unsigned int n0_,
      unsigned int n1_,
      unsigned int n2_,
      unsigned int n3_,
      unsigned int n4_)
      : ptr(ptr_), s0(n1_ * n2_ * n3_ * n4_), s1(n2_ * n3_ * n4_), s2(n3_ * n4_), s3(n4_) {}

  __device__
  T &operator()(
      unsigned int i0,
      unsigned int i1,
      unsigned int i2,
      unsigned int i3,
      unsigned int i4) {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3 * s3 + i4];
  }

  __device__
  const T &operator()(
      unsigned int i0,
      unsigned int i1,
      unsigned int i2,
      unsigned int i3,
      unsigned int i4) const {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3 * s3 + i4];
  }

public:
  const unsigned int s0 = 0, s1 = 0, s2 = 0, s3 = 0;
  T *const ptr;
};

template<typename T>
class ViewTensor5Const {
public:
  __device__
  ViewTensor5Const(const T *ptr_,
                   unsigned int n0_,
                   unsigned int n1_,
                   unsigned int n2_,
                   unsigned int n3_,
                   unsigned int n4_)
      : ptr(ptr_), s0(n1_ * n2_ * n3_ * n4_), s1(n2_ * n3_ * n4_), s2(n3_ * n4_), s3(n4_) {}

  __device__
  const T &operator()(
      unsigned int i0,
      unsigned int i1,
      unsigned int i2,
      unsigned int i3,
      unsigned int i4) const {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3 * s3 + i4];
  }

public:
  const unsigned int s0 = 0, s1 = 0, s2 = 0, s3 = 0;
  const T *const ptr;
};

}
}

#endif