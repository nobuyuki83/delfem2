/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @details OpenGL dependency (not GLUT)
 * this file should not depend anything other than OpenGL
 */

#ifndef DFM2_COLOR_GLOLD_H
#define DFM2_COLOR_GLOLD_H

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cassert>

#include "delfem2/color.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2 {
namespace opengl {

// -------------------------------------------------------------

DFM2_INLINE void myGlMaterialDiffuse(const CColor &color);
DFM2_INLINE void myGlColor(const CColor &c);
DFM2_INLINE void myGlColorDiffuse(const CColor &color);
DFM2_INLINE void myGlDiffuse(const CColor &color);

DFM2_INLINE void DrawBackground(const CColor &c);
DFM2_INLINE void DrawBackground();

DFM2_INLINE void heatmap_glColor(double input);
DFM2_INLINE void heatmap_glDiffuse(double input);

// ------------------------------------------------------------
// std::vector from here

DFM2_INLINE void heatmap(
    double input,
    const std::vector<std::pair<double, CColor> > &colorMap);

DFM2_INLINE void DrawMeshTri3D_VtxColor(
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri,
    std::vector<CColor> &aColor);

// 0: no, 1:lighting, 2:no-lighting
/**
 *
 * @param[in] aXYZ coordinates of the vertex
 * @param[in] aTri triangle index
 * @param[in] aIndGroup store index of group for each triangle. size: "number of triangle"
 * @param[in] aColor
 */
template<typename REAL>
DFM2_INLINE void DrawMeshTri3DFlag_FaceNorm(
    const std::vector<REAL> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<unsigned int> &aFlgElm,
    const std::vector<std::pair<int, CColor> > &aColor);

DFM2_INLINE void DrawMeshTri_ScalarP0(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap);

DFM2_INLINE void DrawMeshTri2D_ScalarP1(
    const double *aXY,
    size_t nXY,
    const unsigned int *aTri,
    size_t nTri,
    const double *paVal,
    int nstride,
    const std::vector<std::pair<double, CColor> > &colorMap);

DFM2_INLINE void DrawMeshTri2D_ScalarP0(
    std::vector<int> &aTri,
    std::vector<double> &aXY,
    std::vector<double> &aVal,
    int nstride,
    int noffset,
    const std::vector<std::pair<double, CColor> > &colorMap);

// 3D value -- vtx value
DFM2_INLINE void DrawMeshTri3D_ScalarP1(
    const double *aXYZ,
    unsigned int nXYZ,
    const unsigned int *aTri,
    unsigned int nTri,
    const double *aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap);

DFM2_INLINE void DrawMeshTri3D_ScalarP1(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const double *aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap);

template<typename REAL>
DFM2_INLINE void DrawMeshTri3D_VtxColor(
    const std::vector<REAL> &aXYZ,
    const std::vector<unsigned int> &aTri,
    std::vector<CColor> &aColor);

DFM2_INLINE void drawMeshTri3D_ScalarP0(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap);

// scalar value on 3D mesh (mixed elem).
DFM2_INLINE void DrawMeshElem3D_Scalar_Vtx(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    const double *aValVtx,
    const std::vector<std::pair<double, CColor> > &colorMap);

// 3D value
DFM2_INLINE void DrawMeshTet3D_ScalarP1(
    const double *aXYZ,
    size_t nXYZ,
    const unsigned int *aTet,
    size_t nTet,
    const double *aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap);

DFM2_INLINE void DrawMeshTet3D_Cut(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTet,
    const std::vector<CColor> &aColor,
    const double org[3], const double n[3]);

} // namespace opengl
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/old/color.cpp"
#endif

#endif
