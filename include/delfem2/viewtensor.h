
#ifndef DFM2_VIEWTENSOR_H
#define DFM2_VIEWTENSOR_H

namespace delfem2 {

template<typename T>
class ViewTensor3 {
public:
  ViewTensor3(
      T *ptr_,
      unsigned int n0_, unsigned int n1_, unsigned int n2_)
      : ptr(ptr_), s0(n1_ * n2_), s1(n2_) {}

  T &operator()(
      unsigned int i0, unsigned int i1, unsigned int i2) {
    return ptr[i0 * s0 + i1 * s1 + i2];
  }

  const T &operator()(
      unsigned int i0, unsigned int i1, unsigned int i2) const {
    return ptr[i0 * s0 + i1 * s1 + i2];
  }

  const T *data(unsigned int i0, unsigned int i1) const {
    return ptr + i0 * s0 + i1 * s1;
  }

public:
  T *const ptr;
  const unsigned int s0 = 0, s1 = 0;
};

template<typename T>
class ViewTensor4 {
public:
  ViewTensor4(
      T *ptr_,
      unsigned int n0_, unsigned int n1_, unsigned int n2_, unsigned int n3_)
      : ptr(ptr_), s0(n1_ * n2_ * n3_), s1(n2_ * n3_), s2(n3_) {}

  T &operator()(
      unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3) {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3];
  }

  const T &operator()(
      unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3) const {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3];
  }

  const T *data(unsigned int i0, unsigned int i1, unsigned int i2) const {
    return ptr + i0 * s0 + i1 * s1 + i2 * s2;
  }

public:
  T *const ptr;
  const unsigned int s0 = 0, s1 = 0, s2 = 0;
};

template<typename T>
class ViewTensor4Const {
public:
  ViewTensor4Const(
      const T *ptr_,
      unsigned int n0_, unsigned int n1_, unsigned int n2_, unsigned int n3_)
      : ptr(ptr_), s0(n1_ * n2_ * n3_), s1(n2_ * n3_), s2(n3_) {}

  const T &operator()(
      unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3) const {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3];
  }

  const T *data(unsigned int i0, unsigned int i1, unsigned int i2) const {
    return ptr + i0 * s0 + i1 * s1 + i2 * s2;
  }

public:
  const T *const ptr;
  const unsigned int s0 = 0, s1 = 0, s2 = 0;
};

template<typename T>
class ViewTensor5Const {
public:
  ViewTensor5Const(
      const T *ptr_,
      unsigned int n0_, unsigned int n1_, unsigned int n2_, unsigned int n3_, unsigned int n4_)
      : ptr(ptr_), s0(n1_ * n2_ * n3_ * n4_), s1(n2_ * n3_ * n4_), s2(n3_ * n4_), s3(n4_) {}

  const T &operator()(
      unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) const {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3 * s3 + i4];
  }

public:
  const T *const ptr;
  const unsigned int s0 = 0, s1 = 0, s2 = 0, s3 = 0;
};

template<typename T>
class ViewTensor5 {
public:
  ViewTensor5(
      T *ptr_,
      unsigned int n0_, unsigned int n1_, unsigned int n2_, unsigned int n3_, unsigned int n4_)
      : ptr(ptr_), s0(n1_ * n2_ * n3_ * n4_), s1(n2_ * n3_ * n4_), s2(n3_ * n4_), s3(n4_) {}

  T &operator()(
      unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3 * s3 + i4];
  }

  const T &operator()(
      unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) const {
    return ptr[i0 * s0 + i1 * s1 + i2 * s2 + i3 * s3 + i4];
  }

public:
  T *const ptr;
  const unsigned int s0 = 0, s1 = 0, s2 = 0, s3 = 0;
};

}


#endif