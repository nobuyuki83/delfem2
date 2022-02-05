//
// Created by Nobuyuki Umetani on 2022/02/04.
//

#ifndef GEO_META_FUNCS_H_
#define GEO_META_FUNCS_H_

namespace delfem2 {

template<class T>
class has_definition {
  template<class VEC>  // VEC is Eigen::Vector3& or dfm2::CVec3&
  static constexpr int check(typename std::remove_reference<VEC>::type::Scalar *) {
    if constexpr( std::is_reference<VEC>::value ){
      return 4;
    }
    return 0;
  }

  template<class VEC>
  static constexpr int check(typename VEC::value_type *) { // VEC is std::array
    static_assert(
        std::tuple_size<VEC>::value == 3,
        "the size of the std::array needs to be 3");
    return 1;
  }
    
  template<class VEC>
  static constexpr int check(typename std::remove_reference_t<VEC>::value_type *) { // VEC is std::array
    using type0 = std::remove_reference_t<VEC>;
    static_assert(
        std::tuple_size<type0>::value == 3,
        "the size of the std::array needs to be 3");
    return 6;
  }

  template<class VEC>
  static constexpr int check(...) {
    if constexpr(std::is_array_v<VEC>) { // this is static array
      return 2;
    }
    else{
        if constexpr(std::is_reference_v<VEC>){ // &&
            return 5;
        }
        else{ // pointer
            static_assert(
          std::is_same_v<std::remove_pointer_t<VEC>*, VEC>,
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

template <class VEC>
struct conditional<4, VEC> {
  using type = typename std::remove_reference_t<VEC>::Scalar;
};

template <class VEC>
struct conditional<5, VEC> {
  using type0 = typename std::remove_reference_t<VEC>;
  using type1 = std::remove_extent_t<type0>;
  using type = std::remove_const_t<type1>;
};

template <class VEC>
struct conditional<6, VEC> {
  using type0 = typename std::remove_reference_t<VEC>;
  using type = typename type0::value_type;
};


template<class T>
using value_type = typename conditional<has_definition<T>::value, T>::type;


}

#endif //GEO_META_FUNCS_H_
