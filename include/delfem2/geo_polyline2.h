/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GEO_POLYLINE2_H
#define DFM2_GEO_POLYLINE2_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

#include "delfem2/vec2.h"
#include "delfem2/dfm2_inline.h"

// -----------------------------------------------------

namespace delfem2 {

// move to paramgeo2d?
template<typename T>
CVec2<T> pointCurve_BezierCubic(
    double t,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    const CVec2<T> &p3,
    const CVec2<T> &p4);

template<typename T>
CVec2<T> pointCurve_BezierQuadratic(
    double t,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    const CVec2<T> &p3);

//! @brief translate all the points
template<typename T>
void Translate(
    std::vector<CVec2<T> > &aP,
    double dx,
    double dy);

template<typename T>
DFM2_INLINE void Rotate(
    std::vector<CVec2<T> > &aP,
    double dt);

template<typename T>
void Polyline_CubicBezierCurve(
    std::vector<CVec2<T> > &aP,
    const int n,
    const std::vector<CVec2<T> > &aCP);

template<typename T>
void Polyline_BezierCubic(
    std::vector<CVec2<T> > &aP,
    const unsigned int n,
    const CVec2<T> &p1, const CVec2<T> &p2,
    const CVec2<T> &p3, const CVec2<T> &p4);

template<typename T>
void Polyline_BezierQuadratic(
    std::vector<CVec2<T>> &aP,
    const unsigned int n,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    const CVec2<T> &p3);

template<typename T>
std::vector<CVec2<T> > Polyline_Resample_Polyline(
    const std::vector<CVec2<T> > &stroke0,
    double l);


} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_polyline2.cpp"
#endif

#endif // VEC_2


