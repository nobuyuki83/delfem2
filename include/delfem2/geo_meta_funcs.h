//
// Created by Nobuyuki Umetani on 2022/02/04.
//

#ifndef GEO_META_FUNCS_H_
#define GEO_META_FUNCS_H_

namespace delfem2 {

template<class T, int N>
class has_definition {

  template<class VEC>
  static constexpr int check(typename VEC::value_type *) { // std::array, Eigen::Vector, CVec3
    return 1;
  }

  template<class VEC>
  static constexpr int check(...) {
    if constexpr(std::is_array_v<VEC>) { // this is static array (e.g., double [3])
      static_assert(
          std::rank_v<VEC> == 1,
          "the dimension of the static array should be one");
      static_assert(
          std::extent_v<VEC,0> == N,
          "the size of the static array is wrong");
      return 2;
    } else {
      static_assert(
          std::is_same_v<std::remove_pointer_t<VEC> *, VEC>,
          "pointer");
      return 3;
    }
    return -1;
  }
 public:
  static constexpr int value = check<T>(nullptr);
};

template<int F, class VEC>
struct conditional {  // VEC is
};

template<class VEC>
struct conditional<1, VEC> {  // VEC is std::array, dfm2::CVec3 or Eigen::Vector3
  using type = typename VEC::value_type;
};

template<class VEC>
struct conditional<2, VEC> {  // VEC is static array (e.g., double [3])
  using type = std::remove_extent_t<VEC>;
};

template<class VEC>
struct conditional<3, VEC> { // VEC is pointer (e.g., double *)
  using type0 = std::remove_pointer_t<VEC>;
  using type = std::remove_const_t<type0>;
};

template<class T, int N, typename T0=std::remove_reference_t<T>>
using vecn_value_t = typename conditional<has_definition<T0,N>::value, T0>::type;

}

#endif //GEO_META_FUNCS_H_
