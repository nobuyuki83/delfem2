/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @details the file depends on OpenGL 2.x.
 * There are a lot of legacy OpenGL commands such as glBegin(),glEnd()
 *
 * (2020/12/26) TODO: split this file into drawedge, drawface.
 */

#ifndef DFM2_OPENGL_OLD_MSHUNI_H
#define DFM2_OPENGL_OLD_MSHUNI_H

#include <string>
#include <vector>

#include "delfem2/dfm2_inline.h"


#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288
#endif

namespace delfem2::opengl {

// ========================================
// Draw Point

DFM2_INLINE void DrawPoints2D_Vectors(
    const double *aXY,
    unsigned int nXY,
    const double *aVal,
    int nstride,
    int noffset,
    double mag);

DFM2_INLINE void DrawPoints2d_Points(
    const std::vector<double> &aXY);

DFM2_INLINE void DrawPoints2d_Psup(
    unsigned int nXY,
    const double *aXY,
    const unsigned int *psup_ind,
    const unsigned int *psup);

DFM2_INLINE void DrawPoints2d_4RotSym(
    const double *aXY,
    unsigned int nXY,
    const double *aDir,
    double vlen);

// above: 2D
// below: 3D

/**
 * @brief Draw Points using GL_POINTS
 * @tparam T float or double
 * @param aXYZ
 */
template<typename T>
DFM2_INLINE void DrawPoints3_Points(
    const std::vector<T> &aXYZ);

DFM2_INLINE void DrawPoints3d_NormVtx(
    const std::vector<double> &aXYZ,
    const std::vector<double> &aNrm,
    double scale);

DFM2_INLINE void DrawPoints3d_Psup(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

// =====================================
// Draw Line

DFM2_INLINE void DrawMeshLine3D_Edge(
    const double *aXYZ,
    unsigned int nXYZ,
    const unsigned int *aLine,
    unsigned int nLine);

// =====================================
// Draw Triangle Mesh

DFM2_INLINE void DrawMeshTri2D_Face(
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXY);

DFM2_INLINE void DrawMeshTri2D_FaceDisp2D(
    const double *aXY,
    size_t nXY,
    const unsigned int *aTri,
    size_t nTri,
    const double *aDisp,
    int nstride);

DFM2_INLINE void DrawMeshTri2D_Edge(
    const double *aXY,
    size_t nXY,
    const unsigned int *aTri,
    size_t nTri);

DFM2_INLINE void DrawMeshTri2D_Edge(
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXY);

DFM2_INLINE void DrawMeshTri2D_FaceColor(
    const unsigned int *aTri,
    size_t nTri,
    const double *aXY,
    size_t nXY,
    const double *aColor);

// above: MeshTri2D
// --------------------
// below: MeshTri3D

DFM2_INLINE void DrawMeshTri3D_FaceEdge(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm(
    const double *paXYZ,
    const unsigned int *paTri,
    size_t nTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aNorm);

DFM2_INLINE void DrawMeshTri3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTriVtx,
    const std::vector<double> &aNorm,
    const std::vector<unsigned int> &aTriNrm);

//
DFM2_INLINE void DrawMeshTri3DPart_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri,
    const std::vector<int> &aIndTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm_Flg(
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri,
    int iflg,
    const std::vector<int> &aFlgTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm_XYsym(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void DrawMeshTri3D_FaceNormEdge(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm_TexFace(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aTex);

DFM2_INLINE void DrawMeshTri3D_FaceNorm_TexVtx(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aTex);

template<typename REAL>
DFM2_INLINE void DrawMeshTri3D_Edge(
    const REAL *aXYZ,
    size_t nXYZ,
    const unsigned int *aTri,
    size_t nTri);

DFM2_INLINE void DrawMeshTri3D_Edge(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void DrawMeshTriMap3D_Edge(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<int> &map);

// =====================================
// Draw Quad Mesh

DFM2_INLINE void DrawMeshQuad2D_Edge(
    const double *aXY,
    size_t nXY,
    const unsigned int *aQuad,
    size_t nQuad);

DFM2_INLINE void DrawMeshQuad2D_Edge(
    const std::vector<double> &aXY,
    const std::vector<unsigned int> &aQuad);

DFM2_INLINE void DrawMeshQuad2D_EdgeDisp(
    const double *aXY,
    unsigned int nXY,
    const unsigned int *aQuad,
    unsigned int nQuad,
    const double *aDisp,
    double scale);

DFM2_INLINE void DrawMeshTri2D_EdgeDisp(
    const double *aXY,
    unsigned int nXY,
    const unsigned int *aTri,
    unsigned int nTri,
    const double *aDisp,
    double scale);

DFM2_INLINE void DrawMeshQuad3D_Edge(
    const double *aXYZ,
    size_t nXYZ,
    const unsigned int *aQuad,
    size_t nQuad);

DFM2_INLINE void DrawMeshQuad3D_Edge(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aQuad);

DFM2_INLINE void DrawMeshQuad3D_FaceNorm(
    const double *aXYZ,
    const unsigned int *aQuad,
    size_t nQuad);

DFM2_INLINE void DrawMeshQuad3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aQuad);

// -------------------
// Draw Tet

DFM2_INLINE void DrawMeshTet3D_Edge(
    const double *aXYZ,
    unsigned int nXYZ,
    const unsigned int *aTet,
    unsigned int nTet);

DFM2_INLINE void DrawMeshTet3D_EdgeDisp(
    const double *aXYZ,
    const unsigned int *aTet, unsigned int nTet,
    const double *aDisp,
    double s0);

DFM2_INLINE void DrawMeshTet3D_FaceNorm(
    const double *aXYZ,
    const unsigned int *aTet,
    unsigned int nTet);

DFM2_INLINE void DrawMeshTet3DSurface_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTet,
    const std::vector<unsigned int> &aTetFace);

DFM2_INLINE void DrawMeshTet3DSurface_Edge(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTet,
    const std::vector<unsigned int> &aTetFace);

DFM2_INLINE void DrawMeshTet3D_FaceNormDisp(
    const double *aXYZ,
    unsigned int nXYZ,
    const unsigned int *aTet,
    unsigned int nTet,
    const double *aDisp);

// -------------
// Draw Hex

DFM2_INLINE void DrawMeshHex3D_Edge(
    const double *aXYZ,
    unsigned int nXYZ,
    const unsigned int *aHex,
    unsigned int nHex);

DFM2_INLINE void DrawMeshHex3D_EdgeDisp(
    const double *vtx_xyz,
    size_t num_vtx,
    const unsigned int *hex_vtx_idx,
    size_t num_hex,
    const double *displacement_xyz_vtx);

DFM2_INLINE void DrawMeshHex3D_FaceNorm(
    const double *aXYZ,
    const unsigned int *aHex,
    unsigned int nHex);

DFM2_INLINE void DrawHex3D_FaceNormDisp(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aHex,
    const std::vector<double> &aDisp);

DFM2_INLINE void Draw_HexMeshFaceDisp(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aHex,
    const std::vector<double> &aDisp);

// -----------
// Draw Mix

DFM2_INLINE void DrawMeshElem3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem);

DFM2_INLINE void DrawMeshElem3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    const std::vector<double> &aUV);

DFM2_INLINE void DrawMeshElemPart3D_FaceNorm_TexPoEl(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    const std::vector<int> &aIndElem,
    const std::vector<double> &aUV);

} // namespace delfem2::opengl

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/old/mshuni.cpp"
#endif

#endif /* utility_gl_h */
