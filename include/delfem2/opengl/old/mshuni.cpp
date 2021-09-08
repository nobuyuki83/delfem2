/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/mshuni.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <climits>

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__) // Mac
#  define GL_SILENCE_DEPRECATION
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

// -----------------------------------------------------------

namespace delfem2::opengl::old::mshuni {

template<int nDim, typename T>
DFM2_INLINE void myGlVtxPtr(const T *p);

template<>
DFM2_INLINE void myGlVtxPtr<3>(const double *p) { ::glVertex3dv(p); }

template<>
DFM2_INLINE void myGlVtxPtr<3>(const float *p) { ::glVertex3fv(p); }

template<>
DFM2_INLINE void myGlVtxPtr<2>(const double *p) { ::glVertex2dv(p); }

template<>
DFM2_INLINE void myGlVtxPtr<2>(const float *p) { ::glVertex2fv(p); }

template<int nno, int ndim>
DFM2_INLINE void FetchDeformedPosition(
    double aR[nno][ndim],
    const unsigned int *aIP,
    const double *aXYZ,
    const double *aDisp,
    double s0,
    int nstride_disp = ndim) {
  for (int ino = 0; ino < nno; ++ino) {
    const unsigned int ip0 = aIP[ino];
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      aR[ino][idim] = aXYZ[ip0 * ndim + idim] + s0 * aDisp[ip0 * nstride_disp + idim];
    }
  }
}

/*
template<int nno, int ndim>
void FetchPosition(
    double* const aP[nno],
    const unsigned int* aIP,
    const double* aPos)
{
  for(int ino=0;ino;++ino) { aP[ino] = aPos + aIP[ino]*ndim; }
}
 */

constexpr int noelElemFace_Hex[8][4] = { // this numbering is corresponds to VTK_HEX
    {0, 4, 7, 3},  // -x
    {1, 2, 6, 5},  // +x
    {0, 1, 5, 4},  // -y
    {3, 7, 6, 2},  // +y
    {0, 3, 2, 1},  // -z
    {4, 5, 6, 7}   // +z
};

constexpr int noelEdge_Hex[12][2] = {
    {0, 1}, {3, 2}, {4, 5}, {7, 6},
    {0, 3}, {1, 2}, {4, 7}, {5, 6},
    {0, 4}, {1, 5}, {3, 7}, {2, 6}
};

DFM2_INLINE void UnitNormalAreaTri3D(
    double n[3], double &a,
    const double v1[3], const double v2[3], const double v3[3]) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  a = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * 0.5;
  const double invlen = 0.5 / a;
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

DFM2_INLINE void Cross3D(
    double r[3],
    const double v1[3],
    const double v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

DFM2_INLINE double Length3D(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

DFM2_INLINE double SquareLength3D(const double v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

DFM2_INLINE void GetVertical2Vector3D(
    const double vec_n[3],
    double vec_x[3],
    double vec_y[3]) {
  const double vec_s[3] = {0, 1, 0};
  Cross3D(vec_x, vec_s, vec_n);
  const double len = Length3D(vec_x);
  if (len < 1.0e-10) {
    const double vec_t[3] = {1, 0, 0};
    Cross3D(vec_x, vec_t, vec_n);  // z????
    Cross3D(vec_y, vec_n, vec_x);  // x????
  } else {
    const double invlen = 1.0 / len;
    vec_x[0] *= invlen;
    vec_x[1] *= invlen;
    vec_x[2] *= invlen;
    Cross3D(vec_y, vec_n, vec_x);
  }
}

DFM2_INLINE void myGlVertex3d(
    unsigned int ixyz,
    const std::vector<double> &aXYZ) {
  ::glVertex3d(aXYZ[ixyz * 3 + 0], aXYZ[ixyz * 3 + 1], aXYZ[ixyz * 3 + 2]);
}

DFM2_INLINE void myGlVertex2d(
    unsigned int ixy,
    const std::vector<double> &aXY) {
  ::glVertex2d(aXY[ixy * 2 + 0], aXY[ixy * 2 + 1]);
}

DFM2_INLINE void myGlNorm3d(
    unsigned int ixyz,
    const std::vector<double> &aNorm) {
  ::glNormal3d(aNorm[ixyz * 3 + 0], aNorm[ixyz * 3 + 1], aNorm[ixyz * 3 + 2]);
}

DFM2_INLINE void DrawSingleTri3D_FaceNorm(
    const double *aXYZ,
    const unsigned int *aIndXYZ,
    const double *pUV) {
  const unsigned int i0 = aIndXYZ[0];  //assert( i0>=0&&i0<(int)aXYZ.size()/3 );
  const unsigned int i1 = aIndXYZ[1];  //assert( i1>=0&&i1<(int)aXYZ.size()/3 );
  const unsigned int i2 = aIndXYZ[2];  //assert( i2>=0&&i2<(int)aXYZ.size()/3 );
  if (i0 == UINT_MAX) {
    assert(i1 == UINT_MAX);
    assert(i2 == UINT_MAX);
    return;
  }
  const double p0[3] = {aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]};
  const double p1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
  const double p2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
  double un[3], area;
  UnitNormalAreaTri3D(un, area, p0, p1, p2);
  ::glNormal3dv(un);
  if (pUV != nullptr) { ::glTexCoord2d(pUV[0], pUV[1]); }
  ::glVertex3dv(p0);
  if (pUV != nullptr) { ::glTexCoord2d(pUV[2], pUV[3]); }
  ::glVertex3dv(p1);
  if (pUV != nullptr) { ::glTexCoord2d(pUV[4], pUV[5]); }
  ::glVertex3dv(p2);
}

DFM2_INLINE void DrawSingleQuad3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *quad_vtx_idx,
    const double *vtx_uv) {
  const unsigned int i0 = quad_vtx_idx[0];  //assert( i0 >= 0 && i0 < (int)aXYZ.size()/3 );
  const unsigned int i1 = quad_vtx_idx[1];  //assert( i1 >= 0 && i1 < (int)aXYZ.size()/3 );
  const unsigned int i2 = quad_vtx_idx[2];  //assert( i2 >= 0 && i2 < (int)aXYZ.size()/3 );
  const unsigned int i3 = quad_vtx_idx[3];  //assert( i3 >= 0 && i3 < (int)aXYZ.size()/3 );
  if (i0 == UINT_MAX) {
    assert(i1 == UINT_MAX && i2 == UINT_MAX && i3 == UINT_MAX);
    return;
  }
  const double p0[3] = {vtx_xyz[i0 * 3 + 0], vtx_xyz[i0 * 3 + 1], vtx_xyz[i0 * 3 + 2]};
  const double p1[3] = {vtx_xyz[i1 * 3 + 0], vtx_xyz[i1 * 3 + 1], vtx_xyz[i1 * 3 + 2]};
  const double p2[3] = {vtx_xyz[i2 * 3 + 0], vtx_xyz[i2 * 3 + 1], vtx_xyz[i2 * 3 + 2]};
  const double p3[3] = {vtx_xyz[i3 * 3 + 0], vtx_xyz[i3 * 3 + 1], vtx_xyz[i3 * 3 + 2]};
  {
    double un0[3], area;
    UnitNormalAreaTri3D(un0, area, p0, p1, p3);
    if (vtx_uv != nullptr) { ::glTexCoord2d(vtx_uv[0], vtx_uv[1]); }
    ::glNormal3dv(un0);
    ::glVertex3dv(p0);
  }
  {
    double un1[3], area;
    UnitNormalAreaTri3D(un1, area, p0, p1, p2);
    if (vtx_uv != nullptr) { ::glTexCoord2d(vtx_uv[2], vtx_uv[3]); }
    ::glNormal3dv(un1);
    ::glVertex3dv(p1);
  }
  {
    double un2[3], area;
    UnitNormalAreaTri3D(un2, area, p1, p2, p3);
    if (vtx_uv != nullptr) { ::glTexCoord2d(vtx_uv[4], vtx_uv[5]); }
    ::glNormal3dv(un2);
    ::glVertex3dv(p2);
  }
  {
    double un3[3], area;
    UnitNormalAreaTri3D(un3, area, p2, p3, p0);
    if (vtx_uv != nullptr) { ::glTexCoord2d(vtx_uv[6], vtx_uv[7]); }
    ::glNormal3dv(un3);
    ::glVertex3dv(p3);
  }
}

DFM2_INLINE void Draw_SurfaceMeshEdge(
    [[maybe_unused]] unsigned int nXYZ,
    const double *paXYZ,
    unsigned int nTri,
    const unsigned int *paTri) {
  ::glEnableClientState(GL_VERTEX_ARRAY);
  ::glVertexPointer(3, GL_DOUBLE, 0, paXYZ);
  ::glBegin(GL_LINES);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const int i1 = static_cast<int>(paTri[itri * 3 + 0]);
    const int i2 = static_cast<int>(paTri[itri * 3 + 1]);
    const int i3 = static_cast<int>(paTri[itri * 3 + 2]);
    ::glArrayElement(i1);
    ::glArrayElement(i2);
    ::glArrayElement(i2);
    ::glArrayElement(i3);
    ::glArrayElement(i3);
    ::glArrayElement(i1);
  }
  ::glEnd();
  ::glDisableClientState(GL_VERTEX_ARRAY);
}

DFM2_INLINE void Draw_SurfaceMeshFace(
    [[maybe_unused]] unsigned int nXYZ,
    const double *paXYZ,
    unsigned int nTri,
    const unsigned int *paTri) {
  ::glEnableClientState(GL_VERTEX_ARRAY);
  ::glVertexPointer(3, GL_DOUBLE, 0, paXYZ);
  ::glDrawElements(
      GL_TRIANGLES,
      static_cast<int>(nTri * 3),
      GL_UNSIGNED_INT,
      paTri);
  ::glDisableClientState(GL_VERTEX_ARRAY);
  /*
    /////
    ::glColor3d(1,1,1);
    ::glBegin(GL_TRIANGLES);
    for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int i1 = paTri[itri*3+0];
    const unsigned int i2 = paTri[itri*3+1];
    const unsigned int i3 = paTri[itri*3+2];
    ::glVertex3dv(paXYZ+i1*3);
    ::glVertex3dv(paXYZ+i2*3);
    ::glVertex3dv(paXYZ+i3*3);
    }
    ::glEnd();
    ::glColor3d(0,0,0);
    ::glBegin(GL_LINES);
    for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int i1 = paTri[itri*3+0];
    const unsigned int i2 = paTri[itri*3+1];
    const unsigned int i3 = paTri[itri*3+2];
    ::glVertex3dv(paXYZ+i1*3);
    ::glVertex3dv(paXYZ+i2*3);
    ::glVertex3dv(paXYZ+i2*3);
    ::glVertex3dv(paXYZ+i3*3);
    ::glVertex3dv(paXYZ+i3*3);
    ::glVertex3dv(paXYZ+i1*3);
    }
    ::glEnd();
    */
}

DFM2_INLINE void drawLoop2d(
    const std::vector<double> &vex_xy) {
  ::glBegin(GL_LINES);
  const int nvec = static_cast<int>(vex_xy.size() / 2);
  for (int ivec = 0; ivec < nvec; ivec++) {
    int jvec = ivec + 1;
    if (jvec >= nvec) { jvec -= nvec; }
    myGlVertex2d(ivec, vex_xy);
    myGlVertex2d(jvec, vex_xy);
  }
  ::glEnd();
  //
  ::glBegin(GL_POINTS);
  for (int ivec = 0; ivec < nvec; ivec++) {
    myGlVertex2d(ivec, vex_xy);
  }
  ::glEnd();
}

DFM2_INLINE void DrawMeshTriMap3D_Edge(
    const std::vector<double> &vtx_xyz,
    const std::vector<int> &tri_vtx_idx,
    const std::vector<int> &map_vtx) {
  const GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  const size_t nTri = tri_vtx_idx.size() / 3;
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const int j1 = tri_vtx_idx[itri * 3 + 0];
    const int j2 = tri_vtx_idx[itri * 3 + 1];
    const int j3 = tri_vtx_idx[itri * 3 + 2];
    if (j1 == -1) {
      assert(j2 == -1);
      assert(j3 == -1);
      continue;
    }
    const int i1 = map_vtx[j1];
    const int i2 = map_vtx[j2];
    const int i3 = map_vtx[j3];
    assert(i1 >= 0 && i1 < (int) vtx_xyz.size() / 3);
    assert(i2 >= 0 && i2 < (int) vtx_xyz.size() / 3);
    assert(i3 >= 0 && i3 < (int) vtx_xyz.size() / 3);
    double p1[3] = {vtx_xyz[i1 * 3 + 0], vtx_xyz[i1 * 3 + 1], vtx_xyz[i1 * 3 + 2]};
    double p2[3] = {vtx_xyz[i2 * 3 + 0], vtx_xyz[i2 * 3 + 1], vtx_xyz[i2 * 3 + 2]};
    double p3[3] = {vtx_xyz[i3 * 3 + 0], vtx_xyz[i3 * 3 + 1], vtx_xyz[i3 * 3 + 2]};
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
    ::glVertex3dv(p3);
    ::glVertex3dv(p1);
  }
  ::glEnd();

  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

DFM2_INLINE void ShadowMatrix(
    float m[16],
    const float plane[4],
    const float lpos[3]) {
  float dot = plane[0] * lpos[0] + plane[1] * lpos[1] + plane[2] * lpos[2] + plane[3];
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      m[j * 4 + i] = -plane[j] * lpos[i];
      if (i == j) { m[j * 4 + i] += dot; }
    }
  }
}

DFM2_INLINE void DrawMeshTri3DPart_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<int> &aIndTri) {
  ::glBegin(GL_TRIANGLES);
  for (int itri : aIndTri) {
    assert(itri >= 0 && itri < (int) aTri.size() / 3);
    DrawSingleTri3D_FaceNorm(aXYZ.data(), aTri.data() + itri * 3, 0);
  }
  ::glEnd();
}

DFM2_INLINE void DrawMeshTri3D_FaceNorm_Flg(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    int iflg,
    const std::vector<int> &aFlgTri) {
  const int nTri = (int) aTri.size() / 3;
  //  const int nXYZ = (int)aXYZ.size()/3;
  ::glBegin(GL_TRIANGLES);
  for (int itri = 0; itri < nTri; ++itri) {
    const int iflg0 = aFlgTri[itri];
    if (iflg0 != iflg) continue;
    DrawSingleTri3D_FaceNorm(aXYZ.data(), aTri.data() + itri * 3, 0);
  }
  ::glEnd();
}

DFM2_INLINE void DrawMeshTri3D_FaceEdge(
    const std::vector<double> &vtx_xyz,
    const std::vector<int> &tri_vtx_idx) {
  const std::size_t nTri = tri_vtx_idx.size() / 3;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const int i1 = tri_vtx_idx[itri * 3 + 0];
    const int i2 = tri_vtx_idx[itri * 3 + 1];
    const int i3 = tri_vtx_idx[itri * 3 + 2];
    myGlVertex3d(i1, vtx_xyz);
    myGlVertex3d(i2, vtx_xyz);
    myGlVertex3d(i3, vtx_xyz);
  }
  ::glEnd();
  // ------------------------------------
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_LINES);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int i1 = tri_vtx_idx[itri * 3 + 0];
    const unsigned int i2 = tri_vtx_idx[itri * 3 + 1];
    const unsigned int i3 = tri_vtx_idx[itri * 3 + 2];
    myGlVertex3d(i1, vtx_xyz);
    myGlVertex3d(i2, vtx_xyz);
    myGlVertex3d(i2, vtx_xyz);
    myGlVertex3d(i3, vtx_xyz);
    myGlVertex3d(i3, vtx_xyz);
    myGlVertex3d(i1, vtx_xyz);
  }
  ::glEnd();
}

}  // namespace delfem2::opengl::old::mshuni

// ==================================================================
// Points

DFM2_INLINE void delfem2::opengl::DrawPoints2D_Vectors(
    const double *vtx_xy,
    size_t num_vtx,
    const double *vtx_displacement,
    int nstride,
    int noffset,
    double scale) {
  ::glBegin(GL_LINES);
  for (unsigned int ino = 0; ino < num_vtx; ino++) {
    const double vx = vtx_displacement[ino * nstride + noffset + 0] * scale;
    const double vy = vtx_displacement[ino * nstride + noffset + 1] * scale;
    const double p0[2] = {vtx_xy[ino * 2 + 0], vtx_xy[ino * 2 + 1]};
    const double p1[2] = {vtx_xy[ino * 2 + 0] + vx, vtx_xy[ino * 2 + 1] + vy};
    ::glVertex2dv(p0);
    ::glVertex2dv(p1);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawPoints2d_Points(
    const std::vector<double> &vtx_xy) {
  const std::size_t nxys = vtx_xy.size() / 2;
  ::glBegin(GL_POINTS);
  for (unsigned int ino = 0; ino < nxys; ino++) {
    const double p0[2] = {vtx_xy[ino * 2 + 0], vtx_xy[ino * 2 + 1]};
    ::glVertex2dv(p0);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawPoints2d_Psup(
    unsigned int num_vtx,
    const double *vtx_xy,
    const unsigned int *psup_ind,
    const unsigned int *psup) {
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip < num_vtx; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      unsigned int jp = psup[ipsup];
      ::glVertex2dv(vtx_xy + ip * 2);
      ::glVertex2dv(vtx_xy + jp * 2);
    }
  }
  ::glEnd();
}

template<typename T>
DFM2_INLINE void
delfem2::opengl::DrawPoints3_Points(
    const std::vector<T> &vtx_xyz) {
  const unsigned int nxyz = vtx_xyz.size() / 3;
  ::glBegin(GL_POINTS);
  for (unsigned int ino = 0; ino < nxyz; ino++) {
    delfem2::opengl::old::mshuni::myGlVtxPtr<3>(vtx_xyz.data() + ino * 3);
  }
  ::glEnd();
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::DrawPoints3_Points(const std::vector<float> &vtx_xyz);
template void delfem2::opengl::DrawPoints3_Points(const std::vector<double> &vtx_xyz);
#endif

// -----------------------------

DFM2_INLINE void
delfem2::opengl::DrawPoints3d_NormVtx(
    const std::vector<double> &aXYZ,
    const std::vector<double> &aNrm,
    double scale) {
  const std::size_t np = aXYZ.size() / 3;
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip < np; ip++) {
    const double p0[3] = {
        aXYZ[ip * 3 + 0],
        aXYZ[ip * 3 + 1],
        aXYZ[ip * 3 + 2]};
    const double p1[3] = {
        aXYZ[ip * 3 + 0] + scale * aNrm[ip * 3 + 0],
        aXYZ[ip * 3 + 1] + scale * aNrm[ip * 3 + 1],
        aXYZ[ip * 3 + 2] + scale * aNrm[ip * 3 + 2]};
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
  }
  ::glEnd();
}

DFM2_INLINE void
delfem2::opengl::DrawPoints3d_Psup(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  const std::size_t np = aXYZ.size() / 3;
  assert(psup_ind.size() == np + 1);
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      unsigned int jp = psup[ipsup];
      ::glVertex3dv(aXYZ.data() + ip * 3);
      ::glVertex3dv(aXYZ.data() + jp * 3);
    }
  }
  ::glEnd();
}

DFM2_INLINE void
delfem2::opengl::DrawPoints2d_4RotSym(
    const double *vtx_xy,
    const unsigned int num_vtx,
    const double *vtx_direction,
    double scale) {
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip < num_vtx; ++ip) {
    ::glColor3d(1, 0, 0);
    ::glVertex2d(vtx_xy[ip * 2 + 0], vtx_xy[ip * 2 + 1]);
    ::glVertex2d(
        vtx_xy[ip * 2 + 0] + scale * vtx_direction[ip * 2 + 0],
        vtx_xy[ip * 2 + 1] + scale * vtx_direction[ip * 2 + 1]);
    //
    ::glColor3d(0, 0, 1);
    ::glVertex2d(vtx_xy[ip * 2 + 0], vtx_xy[ip * 2 + 1]);
    ::glVertex2d(
        vtx_xy[ip * 2 + 0] - scale * vtx_direction[ip * 2 + 0],
        vtx_xy[ip * 2 + 1] - scale * vtx_direction[ip * 2 + 1]);
    //
    ::glVertex2d(
        vtx_xy[ip * 2 + 0] + scale * vtx_direction[ip * 2 + 1],
        vtx_xy[ip * 2 + 1] - scale * vtx_direction[ip * 2 + 0]);
    ::glVertex2d(
        vtx_xy[ip * 2 + 0] - scale * vtx_direction[ip * 2 + 1],
        vtx_xy[ip * 2 + 1] + scale * vtx_direction[ip * 2 + 0]);
  }
  ::glEnd();
}

// ---------------------------------------------
// line

DFM2_INLINE void delfem2::opengl::DrawMeshLine3D_Edge(
    const double *vtx_xyz,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *line_vtx_idx,
    size_t num_line) {
  for (unsigned int il = 0; il < num_line; il++) {
    const unsigned int i0 = line_vtx_idx[il * 2 + 0];
    const unsigned int i1 = line_vtx_idx[il * 2 + 1];
    const double p0[3] = {vtx_xyz[i0 * 3 + 0], vtx_xyz[i0 * 3 + 1], vtx_xyz[i0 * 3 + 2]};
    const double p1[3] = {vtx_xyz[i1 * 3 + 0], vtx_xyz[i1 * 3 + 1], vtx_xyz[i1 * 3 + 2]};
    //::glColor3d(0, 0, 0);
    ::glBegin(GL_LINES);
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glEnd();
  }
}

// =====================================================
// tri3

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *tri_vtx_idx,
    size_t num_tri) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < num_tri; ++itri) {
    lcl::DrawSingleTri3D_FaceNorm(vtx_xyz, tri_vtx_idx + itri * 3, nullptr);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_FaceNorm(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx) {
  DrawMeshTri3D_FaceNorm(
      vtx_xyz.data(),
      tri_vtx_idx.data(),
      tri_vtx_idx.size() / 3);
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_FaceNorm_TexVtx(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aTex) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  const std::size_t nTri = aTri.size() / 3;
  //
  double uv[6];
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const unsigned int ip0 = aTri[itri * 3 + 0];
    const unsigned int ip1 = aTri[itri * 3 + 1];
    const unsigned int ip2 = aTri[itri * 3 + 2];
    uv[0] = aTex[ip0 * 2 + 0];
    uv[1] = aTex[ip0 * 2 + 1];
    uv[2] = aTex[ip1 * 2 + 0];
    uv[3] = aTex[ip1 * 2 + 1];
    uv[4] = aTex[ip2 * 2 + 0];
    uv[5] = aTex[ip2 * 2 + 1];
    lcl::DrawSingleTri3D_FaceNorm(aXYZ.data(), aTri.data() + itri * 3, uv);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_FaceNorm_TexFace(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx,
    const std::vector<double> &vtx_tex) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  const std::size_t nTri = tri_vtx_idx.size() / 3;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    lcl::DrawSingleTri3D_FaceNorm(
        vtx_xyz.data(),
        tri_vtx_idx.data() + itri * 3,
        vtx_tex.data() + itri * 6);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_FaceNorm_XYsym(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  const std::size_t nTri = tri_vtx_idx.size() / 3;
  //
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const unsigned int i1 = tri_vtx_idx[itri * 3 + 0];
    const unsigned int i2 = tri_vtx_idx[itri * 3 + 1];
    const unsigned int i3 = tri_vtx_idx[itri * 3 + 2];
    assert(i1 < vtx_xyz.size() / 3);
    assert(i2 < vtx_xyz.size() / 3);
    assert(i3 < vtx_xyz.size() / 3);
    const double p1[3] = {vtx_xyz[i1 * 3 + 0], vtx_xyz[i1 * 3 + 1], -vtx_xyz[i1 * 3 + 2]};
    const double p2[3] = {vtx_xyz[i2 * 3 + 0], vtx_xyz[i2 * 3 + 1], -vtx_xyz[i2 * 3 + 2]};
    const double p3[3] = {vtx_xyz[i3 * 3 + 0], vtx_xyz[i3 * 3 + 1], -vtx_xyz[i3 * 3 + 2]};
    double un[3], area;
    lcl::UnitNormalAreaTri3D(un, area, p1, p3, p2);
    ::glNormal3dv(un);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    ::glVertex3dv(p2);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_FaceNormEdge(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  const std::size_t nTri = tri_vtx_idx.size() / 3;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    lcl::DrawSingleTri3D_FaceNorm(
        vtx_xyz.data(),
        tri_vtx_idx.data() + itri * 3,
        nullptr);
  }
  ::glEnd();

  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const unsigned int i1 = tri_vtx_idx[itri * 3 + 0];
    const unsigned int i2 = tri_vtx_idx[itri * 3 + 1];
    const unsigned int i3 = tri_vtx_idx[itri * 3 + 2];
    if (i1 == UINT_MAX) {
      assert(i2 == UINT_MAX);
      assert(i3 == UINT_MAX);
      continue;
    }
    assert(i1 < vtx_xyz.size() / 3);
    assert(i2 < vtx_xyz.size() / 3);
    assert(i3 < vtx_xyz.size() / 3);
    const double p1[3] = {vtx_xyz[i1 * 3 + 0], vtx_xyz[i1 * 3 + 1], vtx_xyz[i1 * 3 + 2]};
    const double p2[3] = {vtx_xyz[i2 * 3 + 0], vtx_xyz[i2 * 3 + 1], vtx_xyz[i2 * 3 + 2]};
    const double p3[3] = {vtx_xyz[i3 * 3 + 0], vtx_xyz[i3 * 3 + 1], vtx_xyz[i3 * 3 + 2]};
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
    ::glVertex3dv(p3);
    ::glVertex3dv(p1);
  }
  ::glEnd();

  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

template<typename REAL>
DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_Edge(
    const REAL *vtx_xyz,
    size_t num_vtx,
    const unsigned int *aTri,
    size_t num_tri) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  for (unsigned int itri = 0; itri < num_tri; ++itri) {
    const unsigned int i1 = aTri[itri * 3 + 0];
    const unsigned int i2 = aTri[itri * 3 + 1];
    const unsigned int i3 = aTri[itri * 3 + 2];
    if (i1 >= num_vtx || i2 >= num_vtx || i3 >= num_vtx) { continue; }
    assert(i1 < num_vtx && i2 < num_vtx && i3 < num_vtx);
    const REAL *p1 = vtx_xyz + i1 * 3;
    const REAL *p2 = vtx_xyz + i2 * 3;
    const REAL *p3 = vtx_xyz + i3 * 3;
    lcl::myGlVtxPtr<3>(p1);
    lcl::myGlVtxPtr<3>(p2);
    lcl::myGlVtxPtr<3>(p2);
    lcl::myGlVtxPtr<3>(p3);
    lcl::myGlVtxPtr<3>(p3);
    lcl::myGlVtxPtr<3>(p1);
  }
  ::glEnd();
  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_Edge(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx) {
  DrawMeshTri3D_Edge(
      vtx_xyz.data(),
      vtx_xyz.size() / 3,
      tri_vtx_idx.data(),
      tri_vtx_idx.size() / 3);
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_FaceNorm(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx_idx,
    const std::vector<double> &vtx_normal) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  const std::size_t nTri = tri_vtx_idx.size() / 3;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int i1 = tri_vtx_idx[itri * 3 + 0];
    const unsigned int i2 = tri_vtx_idx[itri * 3 + 1];
    const unsigned int i3 = tri_vtx_idx[itri * 3 + 2];
    lcl::myGlNorm3d(i1, vtx_normal);
    lcl::myGlVertex3d(i1, vtx_xyz);
    lcl::myGlNorm3d(i2, vtx_normal);
    lcl::myGlVertex3d(i2, vtx_xyz);
    lcl::myGlNorm3d(i3, vtx_normal);
    lcl::myGlVertex3d(i3, vtx_xyz);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTriVtx,
    const std::vector<double> &aNorm,
    const std::vector<unsigned int> &aTriNrm) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  const std::size_t nTri = aTriVtx.size() / 3;
  assert(aTriNrm.size() == nTri * 3);
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int iv1 = aTriVtx[itri * 3 + 0];
    const unsigned int iv2 = aTriVtx[itri * 3 + 1];
    const unsigned int iv3 = aTriVtx[itri * 3 + 2];
    const unsigned int in1 = aTriNrm[itri * 3 + 0];
    const unsigned int in2 = aTriNrm[itri * 3 + 1];
    const unsigned int in3 = aTriNrm[itri * 3 + 2];
    const bool bn1 = in1 * 3 < aNorm.size();
    const bool bn2 = in2 * 3 < aNorm.size();
    const bool bn3 = in3 * 3 < aNorm.size();
    const bool bv1 = iv1 * 3 < aXYZ.size();
    const bool bv2 = iv2 * 3 < aXYZ.size();
    const bool bv3 = iv3 * 3 < aXYZ.size();
    const bool bn123 = bn1 && bn2 && bn3;
    const bool bv123 = bv1 && bv2 && bv3;
    if (bn123 && bv123) {
      lcl::myGlNorm3d(in1, aNorm);
      lcl::myGlVertex3d(iv1, aXYZ);
      lcl::myGlNorm3d(in2, aNorm);
      lcl::myGlVertex3d(iv2, aXYZ);
      lcl::myGlNorm3d(in3, aNorm);
      lcl::myGlVertex3d(iv3, aXYZ);
    } else if (bv123) {
      lcl::myGlVertex3d(iv1, aXYZ);
      lcl::myGlVertex3d(iv2, aXYZ);
      lcl::myGlVertex3d(iv3, aXYZ);
    }
  }
  ::glEnd();
}

// tri3
// ------------------------------------------------
// tri2

DFM2_INLINE void delfem2::opengl::DrawMeshTri2D_Face(
    const std::vector<unsigned int> &tri_vtx_idx,
    const std::vector<double> &vtx_xy) {
  const std::size_t ntri = tri_vtx_idx.size() / 3;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < ntri; itri++) {
    const unsigned int i0 = tri_vtx_idx[itri * 3 + 0];
    const unsigned int i1 = tri_vtx_idx[itri * 3 + 1];
    const unsigned int i2 = tri_vtx_idx[itri * 3 + 2];
    const double p0[2] = {vtx_xy[i0 * 2 + 0], vtx_xy[i0 * 2 + 1]};
    const double p1[2] = {vtx_xy[i1 * 2 + 0], vtx_xy[i1 * 2 + 1]};
    const double p2[2] = {vtx_xy[i2 * 2 + 0], vtx_xy[i2 * 2 + 1]};
    ::glVertex2dv(p0);
    ::glVertex2dv(p1);
    ::glVertex2dv(p2);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri2D_FaceDisp2D(
    const double *vtx_xy,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *tri_vtx_idx,
    size_t num_tri,
    const double *vtx_displacement,
    int nstride_disp) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  ::glColor3d(1, 1, 1);
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < num_tri; itri++) {
    const unsigned int *aIP = tri_vtx_idx + itri * 3;
    double aP[3][2];
    lcl::FetchDeformedPosition<3, 2>(aP, aIP, vtx_xy, vtx_displacement, 1., nstride_disp);
    ::glVertex2dv(aP[0]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[2]);
  }
  ::glEnd();
  // --------------------------------------
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_LINES);
  for (unsigned int itri = 0; itri < num_tri; itri++) {
    const unsigned int *aIP = tri_vtx_idx + itri * 3;
    double aP[3][2];
    lcl::FetchDeformedPosition<3, 2>(aP, aIP, vtx_xy, vtx_displacement, 1., nstride_disp);
    ::glVertex2dv(aP[0]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[2]);
    ::glVertex2dv(aP[2]);
    ::glVertex2dv(aP[0]);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri2D_Edge(
    const double *vtx_xy,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *tri_vtx_idx,
    size_t num_tri) {
  //  const unsigned int nxys = (int)aXY.size()/2;
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_LINES);
  for (unsigned int itri = 0; itri < num_tri; itri++) {
    const unsigned int ino0 = tri_vtx_idx[itri * 3 + 0];
    const unsigned int ino1 = tri_vtx_idx[itri * 3 + 1];
    const unsigned int ino2 = tri_vtx_idx[itri * 3 + 2];
    ::glVertex2d(vtx_xy[ino0 * 2 + 0], vtx_xy[ino0 * 2 + 1]);
    ::glVertex2d(vtx_xy[ino1 * 2 + 0], vtx_xy[ino1 * 2 + 1]);
    ::glVertex2d(vtx_xy[ino1 * 2 + 0], vtx_xy[ino1 * 2 + 1]);
    ::glVertex2d(vtx_xy[ino2 * 2 + 0], vtx_xy[ino2 * 2 + 1]);
    ::glVertex2d(vtx_xy[ino2 * 2 + 0], vtx_xy[ino2 * 2 + 1]);
    ::glVertex2d(vtx_xy[ino0 * 2 + 0], vtx_xy[ino0 * 2 + 1]);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri2D_Edge(
    const std::vector<unsigned int> &tri_vtx_idx,
    const std::vector<double> &vtx_xy) {
  DrawMeshTri2D_Edge(
      vtx_xy.data(),
      vtx_xy.size() / 2,
      tri_vtx_idx.data(),
      tri_vtx_idx.size() / 3);
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri2D_FaceColor(
    const unsigned int *tri_vtx_idx,
    size_t num_tri,
    const double *vtx_xy,
    [[maybe_unused]] size_t num_vtx,
    const double *tri_rgb) {
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < num_tri; itri++) {
    const unsigned int i0 = tri_vtx_idx[itri * 3 + 0];
    const unsigned int i1 = tri_vtx_idx[itri * 3 + 1];
    const unsigned int i2 = tri_vtx_idx[itri * 3 + 2];
    ::glColor3dv(tri_rgb + i0 * 3);
    ::glVertex2dv(vtx_xy + i0 * 2);
    ::glColor3dv(tri_rgb + i1 * 3);
    ::glVertex2dv(vtx_xy + i1 * 2);
    ::glColor3dv(tri_rgb + i2 * 3);
    ::glVertex2dv(vtx_xy + i2 * 2);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri2D_EdgeDisp(
    const double *vtx_xy,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *tri_vtx_idx,
    size_t num_tri,
    const double *vtx_displacement,
    double scale) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  // ---------------------
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (unsigned int iq = 0; iq < num_tri; ++iq) {
    const unsigned int *aIP = tri_vtx_idx + iq * 3;
    assert(aIP[0] < num_vtx && aIP[1] < num_vtx && aIP[2] < num_vtx);
    double aP[3][2];
    lcl::FetchDeformedPosition<3, 2>(aP, aIP, vtx_xy, vtx_displacement, scale);
    ::glVertex2dv(aP[0]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[2]);
    ::glVertex2dv(aP[2]);
    ::glVertex2dv(aP[0]);
  }
  ::glEnd();
  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

// above: tri2
// ===============================================================
// below: quad

DFM2_INLINE void delfem2::opengl::DrawMeshQuad3D_Edge(
    const double *vtx_xyz,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *quad_vtx_idx,
    size_t num_quad) {
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (unsigned int iq = 0; iq < num_quad; ++iq) {
    assert(quad_vtx_idx[iq * 4 + 0] < num_vtx);
    assert(quad_vtx_idx[iq * 4 + 1] < num_vtx);
    assert(quad_vtx_idx[iq * 4 + 2] < num_vtx);
    assert(quad_vtx_idx[iq * 4 + 3] < num_vtx);
    const double *aP[4] = {
        vtx_xyz + quad_vtx_idx[iq * 4 + 0] * 3,
        vtx_xyz + quad_vtx_idx[iq * 4 + 1] * 3,
        vtx_xyz + quad_vtx_idx[iq * 4 + 2] * 3,
        vtx_xyz + quad_vtx_idx[iq * 4 + 3] * 3};
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[0]);
  }
  ::glEnd();
  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

DFM2_INLINE void delfem2::opengl::DrawMeshQuad3D_Edge(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &quad_vtx_idx) {
  DrawMeshQuad3D_Edge(
      vtx_xyz.data(),
      vtx_xyz.size() / 3,
      quad_vtx_idx.data(),
      quad_vtx_idx.size() / 4);
}

DFM2_INLINE void delfem2::opengl::DrawMeshQuad3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *quad_vtx_idx,
    const size_t num_quad) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  ::glBegin(GL_QUADS);
  for (unsigned int iq = 0; iq < num_quad; ++iq) {
    lcl::DrawSingleQuad3D_FaceNorm(vtx_xyz, quad_vtx_idx + iq * 4, nullptr);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshQuad3D_FaceNorm(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &quad_vtx_idx) {
  DrawMeshQuad3D_FaceNorm(
      vtx_xyz.data(),
      quad_vtx_idx.data(),
      quad_vtx_idx.size() / 4);
}

DFM2_INLINE void delfem2::opengl::DrawMeshQuad2D_Edge(
    const double *vtx_xy,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *quad_vtx_idx,
    size_t num_quad) {
  const GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  // ---------------------
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (unsigned int iq = 0; iq < num_quad; ++iq) {
    assert(quad_vtx_idx[iq * 4 + 0] < num_vtx);
    assert(quad_vtx_idx[iq * 4 + 1] < num_vtx);
    assert(quad_vtx_idx[iq * 4 + 2] < num_vtx);
    assert(quad_vtx_idx[iq * 4 + 3] < num_vtx);
    const double *aP[4] = {
        vtx_xy + quad_vtx_idx[iq * 4 + 0] * 2,
        vtx_xy + quad_vtx_idx[iq * 4 + 1] * 2,
        vtx_xy + quad_vtx_idx[iq * 4 + 2] * 2,
        vtx_xy + quad_vtx_idx[iq * 4 + 3] * 2};
    ::glVertex2dv(aP[0]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[2]);
    ::glVertex2dv(aP[2]);
    ::glVertex2dv(aP[3]);
    ::glVertex2dv(aP[3]);
    ::glVertex2dv(aP[0]);
  }
  ::glEnd();
  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

DFM2_INLINE void delfem2::opengl::DrawMeshQuad2D_Edge(
    const std::vector<double> &vtx_xy,
    const std::vector<unsigned int> &quad_vtx_idx) {
  DrawMeshQuad2D_Edge(
      vtx_xy.data(),
      vtx_xy.size() / 2,
      quad_vtx_idx.data(),
      quad_vtx_idx.size() / 4);
}

DFM2_INLINE void delfem2::opengl::DrawMeshQuad2D_EdgeDisp(
    const double *vtx_xy,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *quad_vtx_idx,
    size_t num_quad,
    const double *vtx_displacement,
    double scale) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  const GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  // ---------------------
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (unsigned int iq = 0; iq < num_quad; ++iq) {
    const unsigned int *aIP = quad_vtx_idx + iq * 4;
    assert(aIP[0] < num_vtx && aIP[1] < num_vtx && aIP[2] < num_vtx && aIP[3] < num_vtx);
    double aP[4][2];
    lcl::FetchDeformedPosition<4, 2>(aP, aIP, vtx_xy, vtx_displacement, scale);
    ::glVertex2dv(aP[0]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[1]);
    ::glVertex2dv(aP[2]);
    ::glVertex2dv(aP[2]);
    ::glVertex2dv(aP[3]);
    ::glVertex2dv(aP[3]);
    ::glVertex2dv(aP[0]);
  }
  ::glEnd();
  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

// above: quad
// ----------------------------------------------------------------------------
// below: tet

DFM2_INLINE void delfem2::opengl::DrawMeshTet3DSurface_FaceNorm(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tet_vtx_idx,
    const std::vector<unsigned int> &tet_face) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  const unsigned int noelTetFace[4][3] = {
      {1, 2, 3},
      {0, 3, 2},
      {0, 1, 3},
      {0, 2, 1}};
  //  const int nTri = (int)aTri.size()/3;
  // const int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itf = 0; itf < tet_face.size() / 2; ++itf) {
    unsigned int itet = tet_face[itf * 2 + 0];
    unsigned int iface = tet_face[itf * 2 + 1];
    const unsigned int i1 = tet_vtx_idx[itet * 4 + noelTetFace[iface][0]];
    const unsigned int i2 = tet_vtx_idx[itet * 4 + noelTetFace[iface][1]];
    const unsigned int i3 = tet_vtx_idx[itet * 4 + noelTetFace[iface][2]];
    if (i1 == UINT_MAX) {
      assert(i2 == UINT_MAX);
      assert(i3 == UINT_MAX);
      continue;
    }
    assert(i1 < vtx_xyz.size() / 3);
    assert(i2 < vtx_xyz.size() / 3);
    assert(i3 < vtx_xyz.size() / 3);
    double p1[3] = {vtx_xyz[i1 * 3 + 0], vtx_xyz[i1 * 3 + 1], vtx_xyz[i1 * 3 + 2]};
    double p2[3] = {vtx_xyz[i2 * 3 + 0], vtx_xyz[i2 * 3 + 1], vtx_xyz[i2 * 3 + 2]};
    double p3[3] = {vtx_xyz[i3 * 3 + 0], vtx_xyz[i3 * 3 + 1], vtx_xyz[i3 * 3 + 2]};
    double un[3], area;
    lcl::UnitNormalAreaTri3D(un, area, p1, p2, p3);
    ::glNormal3dv(un);
    lcl::myGlVertex3d(i1, vtx_xyz);
    lcl::myGlVertex3d(i2, vtx_xyz);
    lcl::myGlVertex3d(i3, vtx_xyz);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTet3DSurface_Edge(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tet_vtx_idx,
    const std::vector<unsigned int> &vec_tet_face) {
  const unsigned int noelTetFace[4][3] = {
      {1, 2, 3},
      {0, 3, 2},
      {0, 1, 3},
      {0, 2, 1}};
  const GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (unsigned int itf = 0; itf < vec_tet_face.size() / 2; ++itf) {
    const unsigned int itet = vec_tet_face[itf * 2 + 0];
    const unsigned int iface = vec_tet_face[itf * 2 + 1];
    const unsigned int i1 = tet_vtx_idx[itet * 4 + noelTetFace[iface][0]];
    const unsigned int i2 = tet_vtx_idx[itet * 4 + noelTetFace[iface][1]];
    const unsigned int i3 = tet_vtx_idx[itet * 4 + noelTetFace[iface][2]];
    if (i1 == UINT_MAX) {
      assert(i2 == UINT_MAX && i3 == UINT_MAX);
      continue;
    }
    assert(i1 < vtx_xyz.size() / 3);
    assert(i2 < vtx_xyz.size() / 3);
    assert(i3 < vtx_xyz.size() / 3);
    double p1[3] = {vtx_xyz[i1 * 3 + 0], vtx_xyz[i1 * 3 + 1], vtx_xyz[i1 * 3 + 2]};
    double p2[3] = {vtx_xyz[i2 * 3 + 0], vtx_xyz[i2 * 3 + 1], vtx_xyz[i2 * 3 + 2]};
    double p3[3] = {vtx_xyz[i3 * 3 + 0], vtx_xyz[i3 * 3 + 1], vtx_xyz[i3 * 3 + 2]};
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
    ::glVertex3dv(p3);
    ::glVertex3dv(p1);
  }
  ::glEnd();

  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

DFM2_INLINE void delfem2::opengl::DrawMeshTet3D_Edge(
    const double *vtx_xyz,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *tet_vtx_idx,
    size_t num_tet) {
  for (unsigned int itet = 0; itet < num_tet; itet++) {
    const double *aP[4] = {
        vtx_xyz + tet_vtx_idx[itet * 4 + 0] * 3,
        vtx_xyz + tet_vtx_idx[itet * 4 + 1] * 3,
        vtx_xyz + tet_vtx_idx[itet * 4 + 2] * 3,
        vtx_xyz + tet_vtx_idx[itet * 4 + 3] * 3};
    ::glBegin(GL_LINES);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[3]);
    ::glEnd();
  }
}

DFM2_INLINE void delfem2::opengl::DrawMeshTet3D_EdgeDisp(
    const double *vtx_xyz,
    const unsigned int *tet_vtx_idx,
    size_t num_tet,
    const double *vtx_displacement,
    double scale) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  for (unsigned int itet = 0; itet < num_tet; itet++) {
    double aP[4][3];
    lcl::FetchDeformedPosition<4, 3>(
        aP,
        tet_vtx_idx + itet * 4, vtx_xyz, vtx_displacement, scale);
    ::glBegin(GL_LINES);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[3]);
    ::glEnd();
  }
}

DFM2_INLINE void delfem2::opengl::DrawMeshTet3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *tet_vtx_idx,
    size_t num_tet) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  for (unsigned itet = 0; itet < num_tet; itet++) {
    const double *aP[4] = {
        vtx_xyz + tet_vtx_idx[itet * 4 + 0] * 3,
        vtx_xyz + tet_vtx_idx[itet * 4 + 1] * 3,
        vtx_xyz + tet_vtx_idx[itet * 4 + 2] * 3,
        vtx_xyz + tet_vtx_idx[itet * 4 + 3] * 3};
    double un0[3], a0;
    lcl::UnitNormalAreaTri3D(un0, a0, aP[1], aP[2], aP[3]);
    double un1[3], a1;
    lcl::UnitNormalAreaTri3D(un1, a1, aP[2], aP[0], aP[3]);
    double un2[3], a2;
    lcl::UnitNormalAreaTri3D(un2, a2, aP[3], aP[0], aP[1]);
    double un3[3], a3;
    lcl::UnitNormalAreaTri3D(un3, a3, aP[0], aP[2], aP[1]);
    //    ::glColor3d(0, 0, 0);
    ::glBegin(GL_TRIANGLES);
    ::glNormal3dv(un0);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[3]);
    ::glNormal3dv(un1);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[0]);
    ::glNormal3dv(un2);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[1]);
    ::glNormal3dv(un3);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[2]);
    ::glEnd();
  }
}

DFM2_INLINE void delfem2::opengl::DrawMeshTet3D_FaceNormDisp(
    const double *vtx_xyz,
    [[maybe_unused]] const size_t num_vtx,
    const unsigned int *tet_vtx_idx,
    const size_t num_tet,
    const double *vtx_displacement) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  for (unsigned int itet = 0; itet < num_tet; itet++) {
    const unsigned int *aIP = tet_vtx_idx + itet * 4;
    double aP[4][3];
    lcl::FetchDeformedPosition<4, 3>(aP, aIP, vtx_xyz, vtx_displacement, 1);
    double un0[3], a0;
    lcl::UnitNormalAreaTri3D(un0, a0, aP[1], aP[2], aP[3]);
    double un1[3], a1;
    lcl::UnitNormalAreaTri3D(un1, a1, aP[2], aP[0], aP[3]);
    double un2[3], a2;
    lcl::UnitNormalAreaTri3D(un2, a2, aP[3], aP[0], aP[1]);
    double un3[3], a3;
    lcl::UnitNormalAreaTri3D(un3, a3, aP[0], aP[2], aP[1]);
    ::glBegin(GL_TRIANGLES);
    ::glNormal3dv(un0);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[3]);
    ::glNormal3dv(un1);
    ::glVertex3dv(aP[2]);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[0]);
    ::glNormal3dv(un2);
    ::glVertex3dv(aP[3]);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[1]);
    ::glNormal3dv(un3);
    ::glVertex3dv(aP[0]);
    ::glVertex3dv(aP[1]);
    ::glVertex3dv(aP[2]);
    ::glEnd();
  }
}

/*
void opengl::DrawMeshTet3D_FaceNormDisp
 (const double* aXYZ, int nXYZ,
  const unsigned int* aTet, int nTet,
  const double* aDisp)
{
  for (int itet = 0; itet<nTet; itet++){
    const int i0 = aTet[itet*4+0];
    const int i1 = aTet[itet*4+1];
    const int i2 = aTet[itet*4+2];
    const int i3 = aTet[itet*4+3];
    const double p0[3] = { aXYZ[i0*3+0]+aDisp[i0*3+0], aXYZ[i0*3+1]+aDisp[i0*3+1], aXYZ[i0*3+2]+aDisp[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0]+aDisp[i1*3+0], aXYZ[i1*3+1]+aDisp[i1*3+1], aXYZ[i1*3+2]+aDisp[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0]+aDisp[i2*3+0], aXYZ[i2*3+1]+aDisp[i2*3+1], aXYZ[i2*3+2]+aDisp[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0]+aDisp[i3*3+0], aXYZ[i3*3+1]+aDisp[i3*3+1], aXYZ[i3*3+2]+aDisp[i3*3+2] };
    double un0[3], a0; UnitNormalAreaTri3D(un0,a0, p1,p2,p3);
    double un1[3], a1; UnitNormalAreaTri3D(un1,a1, p2,p0,p3);
    double un2[3], a2; UnitNormalAreaTri3D(un2,a2, p3,p0,p1);
    double un3[3], a3; UnitNormalAreaTri3D(un3,a3, p0,p2,p1);
      //    ::glColor3d(0, 0, 0);
    ::glBegin(GL_TRIANGLES);
    ::glNormal3dv(un0); ::glVertex3dv(p1); ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glNormal3dv(un1); ::glVertex3dv(p2); ::glVertex3dv(p3); ::glVertex3dv(p0);
    ::glNormal3dv(un2); ::glVertex3dv(p3); ::glVertex3dv(p0); ::glVertex3dv(p1);
    ::glNormal3dv(un3); ::glVertex3dv(p0); ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glEnd();
  }
}
 */

// above: tet
// ---------------------------------------------------------
// below: hex

DFM2_INLINE void delfem2::opengl::DrawMeshHex3D_FaceNorm(
    const double *vtx_xyz,
    const unsigned int *hex_vtx_idx,
    size_t num_hex) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int ihex = 0; ihex < num_hex; ihex++) {
    const double *aP[8] = {
        vtx_xyz + hex_vtx_idx[ihex * 8 + 0] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 1] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 2] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 3] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 4] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 5] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 6] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 7] * 3};
    for (int iface = 0; iface < 6; ++iface) {
      const double *q0 = aP[lcl::noelElemFace_Hex[iface][0]];
      const double *q1 = aP[lcl::noelElemFace_Hex[iface][1]];
      const double *q2 = aP[lcl::noelElemFace_Hex[iface][2]];
      const double *q3 = aP[lcl::noelElemFace_Hex[iface][3]];
      double un0[3], a0;
      lcl::UnitNormalAreaTri3D(un0, a0, q0, q1, q2);
      ::glNormal3dv(un0);
      ::glVertex3dv(q0);
      ::glVertex3dv(q1);
      ::glVertex3dv(q2);
      double un1[3], a1;
      lcl::UnitNormalAreaTri3D(un1, a1, q0, q2, q3);
      ::glNormal3dv(un1);
      ::glVertex3dv(q0);
      ::glVertex3dv(q2);
      ::glVertex3dv(q3);
    }
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawHex3D_FaceNormDisp(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &hex_vtx_idx,
    const std::vector<double> &vtx_displacement) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int ihex = 0; ihex < hex_vtx_idx.size() / 8; ihex++) {
    const unsigned int *aIP = hex_vtx_idx.data() + ihex * 8;
    double aR[8][3];
    lcl::FetchDeformedPosition<8, 3>(aR, aIP, vtx_xyz.data(), vtx_displacement.data(), 1);
    for (int iface = 0; iface < 6; ++iface) {
      const double *q0 = aR[lcl::noelElemFace_Hex[iface][0]];
      const double *q1 = aR[lcl::noelElemFace_Hex[iface][1]];
      const double *q2 = aR[lcl::noelElemFace_Hex[iface][2]];
      const double *q3 = aR[lcl::noelElemFace_Hex[iface][3]];
      double un0[3], a0;
      lcl::UnitNormalAreaTri3D(un0, a0, q0, q1, q2);
      ::glNormal3dv(un0);
      ::glVertex3dv(q0);
      ::glVertex3dv(q1);
      ::glVertex3dv(q2);
      double un1[3], a1;
      lcl::UnitNormalAreaTri3D(un1, a1, q0, q2, q3);
      ::glNormal3dv(un1);
      ::glVertex3dv(q0);
      ::glVertex3dv(q2);
      ::glVertex3dv(q3);
    }
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshHex3D_Edge(
    const double *vtx_xyz,
    [[maybe_unused]] const size_t num_vtx,
    const unsigned int *hex_vtx_idx,
    const size_t num_hex) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  ::glBegin(GL_LINES);
  for (unsigned int ihex = 0; ihex < num_hex; ihex++) {
    const double *aP[8] = {
        vtx_xyz + hex_vtx_idx[ihex * 8 + 0] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 1] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 2] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 3] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 4] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 5] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 6] * 3,
        vtx_xyz + hex_vtx_idx[ihex * 8 + 7] * 3};
    for (auto iedge : lcl::noelEdge_Hex) {
      ::glVertex3dv(aP[iedge[0]]);
      ::glVertex3dv(aP[iedge[1]]);
    }
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshHex3D_EdgeDisp(
    const double *vtx_xyz,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *hex_vtx_idx,
    size_t num_hex,
    const double *displacement_xyz_vtx) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  ::glBegin(GL_LINES);
  for (unsigned int ihex = 0; ihex < num_hex; ihex++) {
    double aP[8][3];
    lcl::FetchDeformedPosition<8, 3>(
        aP, hex_vtx_idx + ihex * 8,
        vtx_xyz, displacement_xyz_vtx, 1.);
    for (auto iedge : lcl::noelEdge_Hex) {
      ::glVertex3dv(aP[iedge[0]]);
      ::glVertex3dv(aP[iedge[1]]);
    }
  }
  ::glEnd();
}

// above: hex
// ----------------------------------------------------
// below: mix


DFM2_INLINE void delfem2::opengl::DrawMeshElem3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  assert(!aElemInd.empty());
  const std::size_t nelem = aElemInd.size() - 1;
  for (unsigned int ielem = 0; ielem < nelem; ++ielem) {
    const unsigned int ielemind0 = aElemInd[ielem];
    const unsigned int ielemind1 = aElemInd[ielem + 1];
    if (ielemind1 - ielemind0 == 3) {
      ::glBegin(GL_TRIANGLES);
      lcl::DrawSingleTri3D_FaceNorm(
          aXYZ.data(),
          aElem.data() + ielemind0,
          nullptr);
      ::glEnd();
    } else if (ielemind1 - ielemind0 == 4) {
      ::glBegin(GL_QUADS);
      lcl::DrawSingleQuad3D_FaceNorm(
          aXYZ.data(),
          aElem.data() + ielemind0,
          nullptr);
      ::glEnd();
    }
  }
}

DFM2_INLINE void delfem2::opengl::DrawMeshElem3D_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    const std::vector<double> &aUV) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  assert(!aElemInd.empty());
  const std::size_t nelem = aElemInd.size() - 1;
  for (unsigned int ielem = 0; ielem < nelem; ++ielem) {
    const unsigned int ielemind0 = aElemInd[ielem];
    const unsigned int ielemind1 = aElemInd[ielem + 1];
    if (ielemind1 - ielemind0 == 3) {
      ::glBegin(GL_TRIANGLES);
      lcl::DrawSingleTri3D_FaceNorm(
          aXYZ.data(),
          aElem.data() + ielemind0,
          aUV.data() + ielemind0 * 2);
      ::glEnd();
    } else if (ielemind1 - ielemind0 == 4) {
      ::glBegin(GL_QUADS);
      lcl::DrawSingleQuad3D_FaceNorm(
          aXYZ.data(),
          aElem.data() + ielemind0,
          aUV.data() + ielemind0 * 2);
      ::glEnd();
    }
  }
}

DFM2_INLINE void delfem2::opengl::DrawMeshElemPart3D_FaceNorm_TexPoEl(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    const std::vector<int> &aIndElemPart,
    const std::vector<double> &aUV) {
  namespace lcl = ::delfem2::opengl::old::mshuni;
  const bool isUV = (aUV.size() == aElem.size() * 2);
  for (int ielem : aIndElemPart) {
    const unsigned int ielemind0 = aElemInd[ielem];
    const unsigned int ielemind1 = aElemInd[ielem + 1];
    const double *pUV = isUV ? aUV.data() + ielemind0 * 2 : nullptr;
    if (ielemind1 - ielemind0 == 3) {
      ::glBegin(GL_TRIANGLES);
      lcl::DrawSingleTri3D_FaceNorm(
          aXYZ.data(),
          aElem.data() + ielemind0,
          pUV);
      ::glEnd();
    } else if (ielemind1 - ielemind0 == 4) {
      ::glBegin(GL_QUADS);
      lcl::DrawSingleQuad3D_FaceNorm(
          aXYZ.data(),
          aElem.data() + ielemind0,
          pUV);
      ::glEnd();
    }
  }
}
