#ifndef DFM2_VIEW_ARRAY2D_H
#define DFM2_VIEW_ARRAY2D_H

/**
 * view a raw array into array of array such as std::vector<dfm2::CVecXx> std::vector<Eigen::Vector2x>
 * @tparam REAL
 */
template <typename REAL>
class ViewAsArray2D {
 public:
  ViewAsArray2D(const REAL* p, size_t nvtx, size_t ndim)
  : p_(p), nvtx_(nvtx), ndim_(ndim){}

  [[nodiscard]] size_t size() const { return nvtx_; }

  const REAL* operator[](int ivtx) const {
    return p_ + ivtx*ndim_;
  }
 private:
  const REAL* p_;
  const size_t nvtx_;
  const size_t ndim_;
};

#endif