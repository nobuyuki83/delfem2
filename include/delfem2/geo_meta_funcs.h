//
// Created by Nobuyuki Umetani on 2022/02/04.
//

#ifndef GEO_META_FUNCS_H_
#define GEO_META_FUNCS_H_

namespace delfem2 {

template<class T, int N>
class has_definition {

  // VEC is Eigen::Vector3& or dfm2::CVec3&
  template<class VEC>
  static constexpr int check(typename std::remove_reference_t<VEC>::Scalar *) {
    if constexpr(std::is_reference<VEC>::value) {
      return 4;
    }
    return 0;
  }

  template<class VEC>
  static constexpr int check(typename std::remove_reference_t<VEC>::value_type *) { // std::array
    if constexpr(std::is_reference_v<VEC>) {  // VEC is reference of std::array
      using type0 = std::remove_reference_t<VEC>;
      static_assert(
          std::tuple_size<type0>::value == N,
          "the size of the std::array is wrong"  );
      return 6;
    } else {  // VEC is std::array
      static_assert(
          std::tuple_size<VEC>::value == N,
          "the size of the std::array is wrong");
      return 1;
    }

  }

  template<class VEC>
  static constexpr int check(...) {
    if constexpr(std::is_array_v<VEC>) { // this is static array
      static_assert(
          std::rank_v<VEC> == 1,
          "the dimension of the static array is wrong");
      static_assert(
          std::extent_v<VEC,0> == N,
          "the size of the static array is wrong");
      return 2;
    } else {
      if constexpr(std::is_reference_v<VEC>) { // reference
        using type0 = std::remove_reference_t<VEC>;
        if constexpr(std::is_array_v<type0>) {
          static_assert(
              std::rank_v<type0> == 1,
              "the dimension of the static array is wrong");
          static_assert(
              std::extent_v<type0,0> == N,
              "the size of the static array is wrong");
          return 5;
        } else {
          static_assert(
              std::is_pointer_v<type0>,
              "this must be a pointer");
          return 7;
        }
      } else { // pointer
        static_assert(
            std::is_same_v<std::remove_pointer_t<VEC> *, VEC>,
            "pointer");
        return 3;
      }
    }
    return -1;
  }
 public:
  static constexpr int value = check<T>(nullptr);
};

template<int F, class VEC>
struct conditional {  // VEC is dfm2::CVec3 or Eigen::Vector3
  using type = typename VEC::Scalar;
};

template<class VEC>
struct conditional<1, VEC> {  // VEC is std::array
  using type = typename VEC::value_type;
};

template<class VEC>
struct conditional<2, VEC> {  // VEC is static array
  using type = std::remove_extent_t<VEC>;
};

template<class VEC>
struct conditional<3, VEC> { // VEC is pointer
  using type0 = std::remove_pointer_t<VEC>;
  using type = std::remove_const_t<type0>;
};

template<class VEC>
struct conditional<4, VEC> {
  using type = typename std::remove_reference_t<VEC>::Scalar;
};

template<class VEC>
struct conditional<5, VEC> {
  using type0 = std::remove_reference_t<VEC>;
  using type1 = std::remove_extent_t<type0>;
  using type = std::remove_const_t<type1>;
};

template<class VEC>
struct conditional<6, VEC> {
  using type0 = std::remove_reference_t<VEC>;
  using type = typename type0::value_type;
};

template<class VEC>
struct conditional<7, VEC> {
  using type0 = std::remove_reference_t<VEC>;
  using type = std::remove_pointer_t<type0>;
};


template<class T, int N>
using vecn_value_t = typename conditional<has_definition<T,N>::value, T>::type;

}

#endif //GEO_META_FUNCS_H_
