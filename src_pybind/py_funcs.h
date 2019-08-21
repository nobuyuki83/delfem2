/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef py_funcs_h
#define py_funcs_h

#include <cassert>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>


template <typename T>
bool AssertNumpyArray2D
(const pybind11::array_t<T>& np0,
 int n0, int n1)
{
  if( np0.ndim() != 2 ){ return false; }
  if( n0 >= 0 ){ if( np0.shape()[0] != n0 ){ return false; } }
  if( n1 >= 0 ){ if( np0.shape()[1] != n1 ){ return false; } }
  if( np0.strides()[0] != (int)sizeof(T)*np0.shape()[1] ){ return false; }
  if( np0.strides()[1] != (int)sizeof(T) ){ return false; }
  return true;
}

#endif /* py_funcs_h */
