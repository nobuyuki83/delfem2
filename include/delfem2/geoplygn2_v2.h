/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GEOPLYGN2_V2_H
#define DFM2_GEOPLYGN2_V2_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include "delfem2/vec2.h"
#include "delfem2/dfm2_inline.h"

// -----------------------------------------------------

namespace delfem2 {

DFM2_INLINE void makeSplineLoop(
    const std::vector<double>& aCV,
    std::vector<double>& aVecCurve);
  
DFM2_INLINE void MeanValueCoordinate2D(
    double *aW,
    double px, double py,
    const double *aXY, unsigned int nXY);

template <typename T>
std::vector<CVec2<T> > Polygon_Resample_Polygon(
    const std::vector<CVec2<T> >& stroke0,
    double l);
  
template <typename T>
void SecondMomentOfArea_Polygon(
    CVec2<T>& cg,  double& area,
    CVec2<T>& pa1, double& I1,
    CVec2<T>& pa2, double& I2,
    const std::vector<CVec2<T> >& aVec2D);
  
template <typename T>
double Length_Polygon(
    const std::vector<CVec2<T> >& aP);

template <typename T>
double Area_Polygon(
    const std::vector<CVec2<T> >& aP);

template <typename T>
void MeanValueCoordinate(
    std::vector<double>& aW,
    CVec2<T>& p,
    std::vector<CVec2<T> >& aVtx);

template <typename T>
void makeRandomLoop(
    unsigned int nCV,
    std::vector<double>& aCV);

template <typename T>
DFM2_INLINE void makeSplineLoop(
    const std::vector<double>& aCV,
    std::vector<double>& aVecCurve);

template <typename T>
void FixLoopOrientation(
    std::vector<int>& loopIP,
    const std::vector<int>& loopIP_ind,
    const std::vector<CVec2<T> >& aXY);

template <typename T>
std::vector<CVec2<T> > Polygon_Invert(
    const std::vector<CVec2<T> >& aP);

template <typename T>
std::vector<double> XY_Polygon(
    const std::vector<CVec2<T> >& aP);

template <typename T>
void ResamplingLoop(
    std::vector<int>& loopIP1_ind,
    std::vector<int>& loopIP1,
    std::vector<CVec2<T> >& aXY,
    double max_edge_length);

template <typename T>
void JArray_FromVecVec_XY(
    std::vector<int>& aIndXYs,
    std::vector<int>& loopIP0,
    std::vector<CVec2<T> >& aXY,
    const std::vector< std::vector<double> >& aaXY);

template <typename T>
bool IsInclude_Loop(
    const double co[],
    const int ixy_stt, const int ixy_end,
    const std::vector<CVec2<T> >& aXY);

template <typename T>
bool CheckInputBoundaryForTriangulation (
    const std::vector<int>& loopIP_ind,
    const std::vector<CVec2<T> >& aXY);
    
}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/geoplygn2_v2.cpp"
#endif

#endif // VEC_2


