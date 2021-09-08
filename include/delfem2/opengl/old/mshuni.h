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
    const double *vtx_xy,
    size_t num_vtx,
    const double *vtx_displacement,
    int nstride,
    int noffset,
    double scale);

DFM2_INLINE void DrawPoints2d_Points(
    const std::vector<double> &vtx_xy);

DFM2_INLINE void DrawPoints2d_Psup(
    unsigned int num_vtx,
    const double *vtx_xy,
    const unsigned int *psup_ind,
    const unsigned int *psup);

DFM2_INLINE void DrawPoints2d_4RotSym(
    const double *vtx_xy,
    const unsigned int num_vtx,
    const double *vtx_direction,
    double scale);

// above: 2D
// below: 3D

/**
 * @brief Draw Points using GL_POINTS
 * @tparam T float or double
 * @param vtx_xyz
 */
template<typename T>
DFM2_INLINE void DrawPoints3_Points(
    const std::vector<T> &vtx_xyz);

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
    const double *vtx_xyz,
    size_t num_vtx,
    const unsigned int *line_vtx_idx,
    size_t num_line);

// =====================================
// Draw Triangle Mesh

DFM2_INLINE void DrawMeshTri2D_Face(
    const std::vector<unsigned int> &tri_vtx_idx,
    const std::vector<double> &vtx_xy);

DFM2_INLINE void DrawMeshTri2D_FaceDisp2D(
    const double *vtx_xy,
    size_t num_vtx,
    const unsigned int *tri_vtx_idx,
    size_t num_tri,
    const double *vtx_displacement,
    int nstride_disp);

DFM2_INLINE void DrawMeshTri2D_Edge(
    const double *vtx_xy,
    size_t num_vtx,
    const unsigned int *tri_vtx_idx,
    size_t num_tri);

DFM2_INLINE void DrawMeshTri2D_Edge(
    const std::vector<unsigned int> &tri_vtx_idx,
    const std::vector<double> &vtx_xy);

DFM2_INLINE void DrawMeshTri2D_FaceColor(
    const unsigned int *tri_vtx_idx,
    size_t num_tri,
    const double *vtx_xy,
    size_t num_vtx,
    const double *tri_rgb);

// above: MeshTri2D
// --------------------
// below: MeshTri3D

DFM2_INLINE void DrawMeshTri3D_FaceEdge(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *tri_vtx_idx,
    size_t num_tri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx);

DFM2_INLINE void DrawMeshTri3D_FaceNorm(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx,
    const std::vector<double> &vtx_normal);

DFM2_INLINE void DrawMeshTri3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTriVtx,
    const std::vector<double> &aNorm,
    const std::vector<unsigned int> &aTriNrm);

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
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx);

DFM2_INLINE void DrawMeshTri3D_FaceNormEdge(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx);

DFM2_INLINE void DrawMeshTri3D_FaceNorm_TexFace(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx,
    const std::vector<double> &vtx_tex);

DFM2_INLINE void DrawMeshTri3D_FaceNorm_TexVtx(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aTex);

template<typename REAL>
DFM2_INLINE void DrawMeshTri3D_Edge(
    const REAL *vtx_xyz,
    size_t num_vtx,
    const unsigned int *aTri,
    size_t num_tri);

DFM2_INLINE void DrawMeshTri3D_Edge(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx);

DFM2_INLINE void DrawMeshTriMap3D_Edge(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<int> &map);

// =====================================
// Draw Quad Mesh

DFM2_INLINE void DrawMeshQuad2D_Edge(
    const double *vtx_xy,
    size_t num_vtx,
    const unsigned int *quad_vtx_idx,
    size_t num_quad);

DFM2_INLINE void DrawMeshQuad2D_Edge(
    const std::vector<double> &vtx_xy,
    const std::vector<unsigned int> &quad_vtx_idx);

DFM2_INLINE void DrawMeshQuad2D_EdgeDisp(
    const double *vtx_xy,
    size_t num_vtx,
    const unsigned int *quad_vtx_idx,
    size_t num_quad,
    const double *vtx_displacement,
    double scale);

DFM2_INLINE void DrawMeshTri2D_EdgeDisp(
    const double *vtx_xy,
    size_t num_vtx,
    const unsigned int *tri_vtx_idx,
    size_t num_tri,
    const double *vtx_displacement,
    double scale);

DFM2_INLINE void DrawMeshQuad3D_Edge(
    const double *vtx_xyz,
    size_t num_vtx,
    const unsigned int *quad_vtx_idx,
    size_t num_quad);

DFM2_INLINE void DrawMeshQuad3D_Edge(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &quad_vtx_idx);

DFM2_INLINE void DrawMeshQuad3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *quad_vtx_idx,
    const size_t num_quad);

DFM2_INLINE void DrawMeshQuad3D_FaceNorm(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &quad_vtx_idx);

// -------------------
// Draw Tet

DFM2_INLINE void DrawMeshTet3D_Edge(
    const double *vtx_xyz,
    size_t num_vtx,
    const unsigned int *tet_vtx_idx,
    size_t num_tet);

DFM2_INLINE void DrawMeshTet3D_EdgeDisp(
    const double *vtx_xyz,
    const unsigned int *tet_vtx_idx,
    size_t num_tet,
    const double *vtx_displacement,
    double scale);

DFM2_INLINE void DrawMeshTet3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *tet_vtx_idx,
    size_t num_tet);

DFM2_INLINE void DrawMeshTet3DSurface_FaceNorm(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tet_vtx_idx,
    const std::vector<unsigned int> &tet_face);

DFM2_INLINE void DrawMeshTet3DSurface_Edge(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tet_vtx_idx,
    const std::vector<unsigned int> &vec_tet_face);

DFM2_INLINE void DrawMeshTet3D_FaceNormDisp(
    const double *vtx_xyz,
    const size_t num_vtx,
    const unsigned int *tet_vtx_idx,
    const size_t num_tet,
    const double *vtx_displacement);

// -------------
// Draw Hex

DFM2_INLINE void DrawMeshHex3D_Edge(
    const double *vtx_xyz,
    const size_t num_vtx,
    const unsigned int *hex_vtx_idx,
    const size_t num_hex);

DFM2_INLINE void DrawMeshHex3D_EdgeDisp(
    const double *vtx_xyz,
    size_t num_vtx,
    const unsigned int *hex_vtx_idx,
    size_t num_hex,
    const double *displacement_xyz_vtx);

DFM2_INLINE void DrawMeshHex3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *hex_vtx_idx,
    size_t num_hex);

DFM2_INLINE void DrawHex3D_FaceNormDisp(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &hex_vtx_idx,
    const std::vector<double> &vtx_displacement);

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

#endif /* DFM2_OPENGL_OLD_MSHUNI_H */
