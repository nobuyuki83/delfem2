#include <vector>

namespace delfem2 {

template<typename T>
class TensorData2 {
public:
  TensorData2() = default;

  void resize(unsigned int n0_, unsigned int n1_) {
    s0 = n1_;
    aData.resize(n0_ * n1_);
  }

  void resize(unsigned int n0_, unsigned int n1_, T v) {
    s0 = n1_;
    aData.resize(n0_ * n1_, v);
  }

  T &operator()(unsigned int i0, unsigned int i1) {
    return aData[i0 * s0 + i1];
  }

  const T &operator()(unsigned int i0, unsigned int i1) const {
    return aData[i0 * s0 + i1];
  }

public:
  unsigned int s0 = 0;
  std::vector <T> aData;
};

template<typename T>
class TensorData3 {
public:
  TensorData3() = default;

  void resize(unsigned int n0_, unsigned int n1_, unsigned int n2_) {
    s0 = n2_ * n1_;
    s1 = n2_;
    aData.resize(n0_ * n1_ * n2_);
  }

  void resize(unsigned int n0_, unsigned int n1_, unsigned int n2_, T v) {
    s0 = n2_ * n1_;
    s1 = n2_;
    aData.resize(n0_ * n1_ * n2_, v);
  }

  const T *data(unsigned int i0, unsigned int i1) const {
    return aData.data() + i0 * s0 + i1 * s1;
  }

  T *data(unsigned int i0, unsigned int i1) {
    return aData.data() + i0 * s0 + i1 * s1;
  }

  T &operator()(unsigned int i0, unsigned int i1, unsigned int i2) {
    return aData[i0 * s0 + i1 * s1 + i2];
  }

  const T &operator()(unsigned int i0, unsigned int i1, unsigned int i2) const {
    return aData[i0 * s0 + i1 * s1 + i2];
  }

public:
  unsigned int s0 = 0, s1 = 0;
  std::vector <T> aData;
};

template<typename T>
class TensorData4 {
public:
  TensorData4() = default;

  void resize(unsigned int n0_, unsigned int n1_, unsigned int n2_, unsigned int n3_) {
    s0 = n3_ * n2_ * n1_;
    s1 = n3_ * n2_;
    s2 = n3_;
    aData.resize(n0_ * n1_ * n2_ * n3_);
  }

  const T *data(unsigned int i0, unsigned int i1, unsigned int i2) const {
    return aData.data() + i0 * s0 + i1 * s1 + i2 * s2;
  }

  T *data(unsigned int i0, unsigned int i1, unsigned int i2) {
    return aData.data() + i0 * s0 + i1 * s1 + i2 * s2;
  }

  T &operator()(unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3) {
    return aData[i0 * s0 + i1 * s1 + i2 * s2 + i3];
  }

  const T &operator()(unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3) const {
    return aData[i0 * s0 + i1 * s1 + i2 * s2 + i3];
  }

public:
  unsigned int s0 = 0, s1 = 0, s2 = 0;
  std::vector <T> aData;
};


template<typename T>
class TensorData5 {
public:
  TensorData5() = default;

  void resize(unsigned int n0_, unsigned int n1_, unsigned int n2_, unsigned int n3_, unsigned int n4_) {
    s0 = n4_ * n3_ * n2_ * n1_;
    s1 = n4_ * n3_ * n2_;
    s2 = n4_ * n3_;
    s3 = n4_;
    aData.resize(n0_ * n1_ * n2_ * n3_ * n4_);
  }

  T &operator()(unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) {
    return aData[i0 * s0 + i1 * s1 + i2 * s2 + i3 * s3 + i4];
  }

  const T &operator()(unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) const {
    return aData[i0 * s0 + i1 * s1 + i2 * s2 + i3 * s3 + i4];
  }

public:
  unsigned int s0 = 0, s1 = 0, s2 = 0, s3 = 0;
  std::vector <T> aData;
};


}