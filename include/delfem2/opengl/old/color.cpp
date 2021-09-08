/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/color.h"

#include <climits>
#include <utility>  // for std::pair

#if defined(_WIN32)  // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
#  define GL_SILENCE_DEPRECATION
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

#include "delfem2/color.h"

#if defined(_MSC_VER)
#  pragma warning( push )
#  pragma warning( disable : 4100 )
#endif

// header ends here
// -------------------------------------------------

namespace delfem2::opengl::color_glold {

template<typename REAL>
DFM2_INLINE void UnitNormalAreaTri3D(
    REAL n[3],
    REAL &a,
    const REAL v1[3],
    const REAL v2[3],
    const REAL v3[3]) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  a = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * 0.5;
  const REAL invlen = 0.5 / a;
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

template<typename T>
void myGlNormalPtr(const T *p);

template<>
void myGlNormalPtr(const double *p) { ::glNormal3dv(p); }

template<>
void myGlNormalPtr(const float *p) { ::glNormal3fv(p); }

template<typename REAL>
DFM2_INLINE void myGlVertex3(
    unsigned int i,
    const std::vector<REAL> &aV);

template<>
DFM2_INLINE void myGlVertex3(
    unsigned int i,
    const std::vector<float> &aV) {
  glVertex3f(aV[i * 3 + 0], aV[i * 3 + 1], aV[i * 3 + 2]);
}

template<>
DFM2_INLINE void myGlVertex3(
    unsigned int i,
    const std::vector<double> &aV) {
  glVertex3d(aV[i * 3 + 0], aV[i * 3 + 1], aV[i * 3 + 2]);
}

DFM2_INLINE void DrawSingleTri3D_Scalar_Vtx(
    const double *aXYZ,
    const unsigned int *tri,
    const double *aValVtx,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  const unsigned int i0 = tri[0];
  const unsigned int i1 = tri[1];
  const unsigned int i2 = tri[2];
  if (i0 == UINT_MAX) {
    assert(i1 == UINT_MAX && i2 == UINT_MAX);
    return;
  }
  //  assert(i0>=0&&i0<(int)aXYZ.size()/3);
  //  assert(i1>=0&&i1<(int)aXYZ.size()/3);
  //  assert(i2>=0&&i2<(int)aXYZ.size()/3);
  const double p0[3] = {aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]};
  const double p1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
  const double p2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
  {
    double n[3], a;
    UnitNormalAreaTri3D(n, a, p0, p1, p2);
    ::glNormal3dv(n);
  }
  const double vt0 = aValVtx[i0];
  const double vt1 = aValVtx[i1];
  const double vt2 = aValVtx[i2];
  heatmap(vt0, colorMap);
  glVertex3dv(p0);
  heatmap(vt1, colorMap);
  glVertex3dv(p1);
  heatmap(vt2, colorMap);
  glVertex3dv(p2);
}

DFM2_INLINE void DrawSingleQuad3D_Scalar_Vtx(
    const std::vector<double> &aXYZ,
    const unsigned int *quad,
    const double *aValVtx,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  const unsigned int i0 = quad[0];
  const unsigned int i1 = quad[1];
  const unsigned int i2 = quad[2];
  const unsigned int i3 = quad[3];
  assert(i0 < aXYZ.size() / 3);
  assert(i1 < aXYZ.size() / 3);
  assert(i2 < aXYZ.size() / 3);
  assert(i3 < aXYZ.size() / 3);
  const double p0[3] = {aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]};
  const double p1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
  const double p2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
  const double p3[3] = {aXYZ[i3 * 3 + 0], aXYZ[i3 * 3 + 1], aXYZ[i3 * 3 + 2]};
  {
    double n[3], a;
    UnitNormalAreaTri3D(n, a, p0, p1, p2);
    ::glNormal3dv(n);
  }
  const double vt0 = aValVtx[i0];
  const double vt1 = aValVtx[i1];
  const double vt2 = aValVtx[i2];
  const double vt3 = aValVtx[i3];
  heatmap(vt0, colorMap);
  glVertex3dv(p0);
  heatmap(vt1, colorMap);
  glVertex3dv(p1);
  heatmap(vt2, colorMap);
  glVertex3dv(p2);
  heatmap(vt3, colorMap);
  glVertex3dv(p3);
}

DFM2_INLINE bool IsAbovePlane(
    const double p[3],
    const double org[3],
    const double n[3]) {
  const double dot =
      (p[0] - org[0]) * n[0] +
          (p[1] - org[1]) * n[1] +
          (p[2] - org[2]) * n[2];
  return dot > 0;
}

} // namespace delfem2::opengl::color_glold

// -----------------------------------------------------------------

DFM2_INLINE void delfem2::opengl::myGlMaterialDiffuse(
    const CColor &color) {
  float c[4];
  c[0] = color.r;
  c[1] = color.g;
  c[2] = color.b;
  c[3] = color.a;
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

DFM2_INLINE void delfem2::opengl::myGlColor(
    const CColor &c) {
  ::glColor4d(c.r, c.g, c.b, c.a);
}

DFM2_INLINE void delfem2::opengl::myGlColorDiffuse(
    const CColor &color) {
  ::glColor4d(color.r, color.g, color.b, color.a);
  float c[4] = {color.r, color.g, color.b, color.a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

DFM2_INLINE void delfem2::opengl::myGlDiffuse(
    const CColor &color) {
  float c[4] = {color.r, color.g, color.b, color.a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

DFM2_INLINE void delfem2::opengl::DrawBackground(
    const CColor &c) {
  glPushAttrib(GL_TRANSFORM_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT);
  ::glShadeModel(GL_SMOOTH);
  GLboolean is_lighting = (glIsEnabled(GL_LIGHTING));
  GLboolean is_texture = (glIsEnabled(GL_TEXTURE_2D));
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_TEXTURE_2D);
  ::glDisable(GL_DEPTH_TEST);

  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glLoadIdentity();
  ::glMatrixMode(GL_PROJECTION);
  ::glPushMatrix();
  ::glLoadIdentity();

  ::glBegin(GL_QUADS);
  ::glColor3f(c.r, c.g, c.b);
  ::glVertex3d(-1, -1, 0);
  ::glVertex3d(+1, -1, 0);
  ::glColor3d(1, 1, 1);
  ::glVertex3d(+1, +1, 0);
  ::glVertex3d(-1, +1, 0);
  ::glEnd();

  ::glMatrixMode(GL_PROJECTION);
  ::glPopMatrix();
  ::glMatrixMode(GL_MODELVIEW);
  ::glPopMatrix();

  ::glPopAttrib();

  ::glEnable(GL_DEPTH_TEST);
  if (is_lighting) { ::glEnable(GL_LIGHTING); }
  if (is_texture) { ::glEnable(GL_TEXTURE_2D); }
}

DFM2_INLINE void delfem2::opengl::DrawBackground() {
  //  ::glColor3d(0.2,0.7,0.7);
  opengl::DrawBackground(CColor(0.5, 0.5, 0.5));
}

// ------------------------------------------------------------

DFM2_INLINE void delfem2::opengl::heatmap_glColor(
    double input) {
  double c[3];
  ::delfem2::heatmap(input, c);
  ::glColor3dv(c);
}

DFM2_INLINE void delfem2::opengl::heatmap_glDiffuse(
    double input) {
  double c[3];
  ::delfem2::heatmap(input, c);
  float cf[4] = {
      (float) c[0],
      (float) c[1],
      (float) c[2],
      1.f};
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, cf);
}

DFM2_INLINE void delfem2::opengl::heatmap(
    double input,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  const CColor &c = getColor(input, colorMap);
  opengl::myGlColorDiffuse(c);
}

// -------------------------------------------------------------

template<typename REAL>
DFM2_INLINE void delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(
    const std::vector<REAL> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<unsigned int> &aFlgElm,
    const std::vector<std::pair<int, CColor> > &aColor) {
  namespace lcl = delfem2::opengl::color_glold;
  const size_t nTri = aTri.size() / 3;
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const unsigned int ig0 = aFlgElm[itri];
    if (ig0 >= aColor.size()) continue;
    const int imode = aColor[ig0].first;
    if (imode == 0) continue;
    else if (imode == 1) { ::glEnable(GL_LIGHTING); }
    else if (imode == 2) { ::glDisable(GL_LIGHTING); }
    myGlColorDiffuse(aColor[ig0].second);
    const unsigned int i1 = aTri[itri * 3 + 0];
    const unsigned int i2 = aTri[itri * 3 + 1];
    const unsigned int i3 = aTri[itri * 3 + 2];
    if (i1 == UINT_MAX) {
      assert(i2 == UINT_MAX && i3 == UINT_MAX);
      continue;
    }
    ::glBegin(GL_TRIANGLES);
    assert(i1 < aXYZ.size() / 3);
    assert(i2 < aXYZ.size() / 3);
    assert(i3 < aXYZ.size() / 3);
    const REAL *p1 = aXYZ.data() + i1 * 3;
    const REAL *p2 = aXYZ.data() + i2 * 3;
    const REAL *p3 = aXYZ.data() + i3 * 3;
    REAL un[3], area;
    lcl::UnitNormalAreaTri3D(un, area, p1, p2, p3);
    lcl::myGlNormalPtr(un);
    lcl::myGlVertex3<REAL>(i1, aXYZ);
    lcl::myGlVertex3<REAL>(i2, aXYZ);
    lcl::myGlVertex3<REAL>(i3, aXYZ);
    ::glEnd();
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(
    const std::vector<float> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<unsigned int> &aFlgElm,
    const std::vector<std::pair<int, CColor> > &aColor);
template void delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<unsigned int> &aFlgElm,
    const std::vector<std::pair<int, CColor> > &aColor);
#endif

DFM2_INLINE void delfem2::opengl::DrawMeshTri2D_ScalarP1(
    const double *aXY,
    size_t nXY,
    const unsigned int *aTri,
    size_t nTri,
    const double *paVal,
    int nstride,
    const std::vector<std::pair<double, CColor> > &colorMap) {
//  const unsigned int ntri = (int)aTri.size()/3;
//  const unsigned int nxys = (int)aXY.size()/2;
  glShadeModel(GL_SMOOTH);
  ::glColor3d(1, 1, 1);
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const unsigned int ino0 = aTri[itri * 3 + 0];
    const unsigned int ino1 = aTri[itri * 3 + 1];
    const unsigned int ino2 = aTri[itri * 3 + 2];
    assert(ino0 < nXY);
    assert(ino1 < nXY);
    assert(ino2 < nXY);
    const double v0 = paVal[ino0 * nstride];
    const double v1 = paVal[ino1 * nstride];
    const double v2 = paVal[ino2 * nstride];
    opengl::heatmap(v0, colorMap);
    ::glVertex2d(aXY[ino0 * 2 + 0], aXY[ino0 * 2 + 1]);
    opengl::heatmap(v1, colorMap);
    ::glVertex2d(aXY[ino1 * 2 + 0], aXY[ino1 * 2 + 1]);
    opengl::heatmap(v2, colorMap);
    ::glVertex2d(aXY[ino2 * 2 + 0], aXY[ino2 * 2 + 1]);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri2D_ScalarP0(
    std::vector<int> &aTri,
    std::vector<double> &aXY,
    std::vector<double> &aVal,
    int nstride,
    int noffset,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  const size_t ntri = aTri.size() / 3;
  ::glColor3d(1, 1, 1);
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < ntri; ++itri) {
    const int ino0 = aTri[itri * 3 + 0];
    const int ino1 = aTri[itri * 3 + 1];
    const int ino2 = aTri[itri * 3 + 2];
    const double v0 = aVal[itri * nstride + noffset];
    opengl::heatmap(v0, colorMap);
    ::glVertex2d(aXY[ino0 * 2 + 0], aXY[ino0 * 2 + 1]);
    ::glVertex2d(aXY[ino1 * 2 + 0], aXY[ino1 * 2 + 1]);
    ::glVertex2d(aXY[ino2 * 2 + 0], aXY[ino2 * 2 + 1]);
  }
  ::glEnd();
}

// vetex value
DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_ScalarP1(
    const double *aXYZ,
    unsigned int nXYZ,
    const unsigned int *aTri,
    unsigned int nTri,
    const double *aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    color_glold::DrawSingleTri3D_Scalar_Vtx(
        aXYZ, aTri + itri * 3, aValSrf, colorMap);
  }
  ::glEnd();
}

// vetex value
DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_ScalarP1(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const double *aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  const size_t nTri = aTri.size() / 3;
  const size_t nXYZ = aXYZ.size() / 3;
  DrawMeshTri3D_ScalarP1(
      aXYZ.data(), static_cast<unsigned int>(nXYZ),
      aTri.data(), static_cast<unsigned int>(nTri),
      aValSrf,
      colorMap);
}

// vetex value
DFM2_INLINE void delfem2::opengl::DrawMeshElem3D_Scalar_Vtx(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    const double *aValVtx,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  if (aElemInd.empty()) return;
  //
  for (size_t ielem = 0; ielem < aElemInd.size() - 1; ++ielem) {
    const unsigned int ielemind0 = aElemInd[ielem];
    const unsigned int ielemind1 = aElemInd[ielem + 1];
    if (ielemind1 - ielemind0 == 3) {
      ::glBegin(GL_TRIANGLES);
      color_glold::DrawSingleTri3D_Scalar_Vtx(
          aXYZ.data(),
          aElem.data() + ielemind0,
          aValVtx,
          colorMap);
      ::glEnd();
    } else if (ielemind1 - ielemind0 == 4) {
      ::glBegin(GL_QUADS);
      color_glold::DrawSingleQuad3D_Scalar_Vtx(
          aXYZ,
          aElem.data() + ielemind0,
          aValVtx,
          colorMap);
      ::glEnd();
    }
  }
}

// element-wise
DFM2_INLINE void delfem2::opengl::drawMeshTri3D_ScalarP0(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  const size_t nTri = aTri.size() / 3;
  if (aValSrf.size() != nTri) { return; }
  // ---
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const unsigned int i0 = aTri[itri * 3 + 0];
    const unsigned int i1 = aTri[itri * 3 + 1];
    const unsigned int i2 = aTri[itri * 3 + 2];
    if (i0 == UINT_MAX) {
      assert(i1 == UINT_MAX);
      assert(i2 == UINT_MAX);
      continue;
    }
    assert(i0 < aXYZ.size() / 3);
    assert(i1 < aXYZ.size() / 3);
    assert(i2 < aXYZ.size() / 3);
    const double p0[3] = {aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]};
    const double p1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
    const double p2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
    double n[3], a;
    color_glold::UnitNormalAreaTri3D(n, a, p0, p1, p2);
    const double vt = aValSrf[itri];
    ::glNormal3dv(n);
    heatmap(vt, colorMap);
    glVertex3dv(p0);
    glVertex3dv(p1);
    glVertex3dv(p2);
  }
  ::glEnd();
}

template<typename REAL>
DFM2_INLINE void delfem2::opengl::DrawMeshTri3D_VtxColor(
    const std::vector<REAL> &aXYZ,
    const std::vector<unsigned int> &aTri,
    std::vector<CColor> &aColor) {
  const unsigned int nTri = aTri.size() / 3;
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const unsigned int i1 = aTri[itri * 3 + 0];
    const unsigned int i2 = aTri[itri * 3 + 1];
    const unsigned int i3 = aTri[itri * 3 + 2];
    if (i1 == UINT_MAX) {
      assert(i2 == UINT_MAX && i3 == UINT_MAX);
      continue;
    }
    ::glBegin(GL_TRIANGLES);
    assert(i1 < aXYZ.size() / 3);
    assert(i2 < aXYZ.size() / 3);
    assert(i3 < aXYZ.size() / 3);
    const REAL p1 = aXYZ.data() + i1 * 3;
    const REAL p2 = aXYZ.data() + i2 * 3;
    const REAL p3 = aXYZ.data() + i3 * 3;
    double un[3], area;
    color_glold::UnitNormalAreaTri3D(un, area, p1, p2, p3);
    ::glNormal3dv(un);
    myGlColorDiffuse(aColor[i1]);
    color_glold::myGlVertex3<double>(i1, aXYZ);
    myGlColorDiffuse(aColor[i2]);
    color_glold::myGlVertex3<double>(i2, aXYZ);
    myGlColorDiffuse(aColor[i3]);
    color_glold::myGlVertex3<double>(i3, aXYZ);
    ::glEnd();
  }
}



// -------------------------------------------------
// tet from here

// 3D value
DFM2_INLINE void delfem2::opengl::DrawMeshTet3D_ScalarP1(
    const double *aXYZ,
    size_t nXYZ,
    const unsigned int *aTet,
    size_t nTet,
    const double *aValSrf,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  ::glBegin(GL_TRIANGLES);
  for (unsigned itri = 0; itri < nTet; ++itri) {
    const unsigned int i0 = aTet[itri * 4 + 0];
    assert(i0 < nXYZ);
    const unsigned int i1 = aTet[itri * 4 + 1];
    assert(i1 < nXYZ);
    const unsigned int i2 = aTet[itri * 4 + 2];
    assert(i2 < nXYZ);
    const unsigned int i3 = aTet[itri * 4 + 3];
    assert(i3 < nXYZ);
    const double p0[3] = {aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]};
    const double p1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
    const double p2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
    const double p3[3] = {aXYZ[i3 * 3 + 0], aXYZ[i3 * 3 + 1], aXYZ[i3 * 3 + 2]};
    double un0[3], a0;
    color_glold::UnitNormalAreaTri3D(un0, a0, p1, p2, p3);
    double un1[3], a1;
    color_glold::UnitNormalAreaTri3D(un1, a1, p2, p0, p3);
    double un2[3], a2;
    color_glold::UnitNormalAreaTri3D(un2, a2, p3, p0, p1);
    double un3[3], a3;
    color_glold::UnitNormalAreaTri3D(un3, a3, p0, p2, p1);
    const double vt0 = aValSrf[i0];
    const double vt1 = aValSrf[i1];
    const double vt2 = aValSrf[i2];
    const double vt3 = aValSrf[i3];
    ::glNormal3dv(un0);
    heatmap(vt1, colorMap);
    glVertex3dv(p1);
    heatmap(vt2, colorMap);
    glVertex3dv(p2);
    heatmap(vt3, colorMap);
    glVertex3dv(p3);
    ::glNormal3dv(un1);
    heatmap(vt2, colorMap);
    glVertex3dv(p2);
    heatmap(vt3, colorMap);
    glVertex3dv(p3);
    heatmap(vt0, colorMap);
    glVertex3dv(p0);
    ::glNormal3dv(un2);
    heatmap(vt3, colorMap);
    glVertex3dv(p3);
    heatmap(vt0, colorMap);
    glVertex3dv(p0);
    heatmap(vt1, colorMap);
    glVertex3dv(p1);
    ::glNormal3dv(un3);
    heatmap(vt0, colorMap);
    glVertex3dv(p0);
    heatmap(vt1, colorMap);
    glVertex3dv(p1);
    heatmap(vt2, colorMap);
    glVertex3dv(p2);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTet3D_Cut(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTet,
    const std::vector<CColor> &aColor,
    const double org[3],
    const double ncut[3]) {
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_DIFFUSE);
  ::glColor3d(1, 1, 1);
  ::glBegin(GL_TRIANGLES);
  for (size_t itet = 0; itet < aTet.size() / 4; itet++) {
    const unsigned int ino0 = aTet[itet * 4 + 0];
    const unsigned int ino1 = aTet[itet * 4 + 1];
    const unsigned int ino2 = aTet[itet * 4 + 2];
    const unsigned int ino3 = aTet[itet * 4 + 3];
    const double p0[3] = {aXYZ[ino0 * 3 + 0], aXYZ[ino0 * 3 + 1], aXYZ[ino0 * 3 + 2]};
    const double p1[3] = {aXYZ[ino1 * 3 + 0], aXYZ[ino1 * 3 + 1], aXYZ[ino1 * 3 + 2]};
    const double p2[3] = {aXYZ[ino2 * 3 + 0], aXYZ[ino2 * 3 + 1], aXYZ[ino2 * 3 + 2]};
    const double p3[3] = {aXYZ[ino3 * 3 + 0], aXYZ[ino3 * 3 + 1], aXYZ[ino3 * 3 + 2]};
    if (color_glold::IsAbovePlane(p0, org, ncut)) continue;
    if (color_glold::IsAbovePlane(p1, org, ncut)) continue;
    if (color_glold::IsAbovePlane(p2, org, ncut)) continue;
    if (color_glold::IsAbovePlane(p3, org, ncut)) continue;
    //    ::glColor3d(1,1,0);
    opengl::myGlColorDiffuse(aColor[itet]);
    //
    double n[3], area;
    color_glold::UnitNormalAreaTri3D(n, area, p0, p2, p1);
    ::glNormal3dv(n);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    //
    color_glold::UnitNormalAreaTri3D(n, area, p0, p1, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    //
    color_glold::UnitNormalAreaTri3D(n, area, p1, p2, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
    //
    color_glold::UnitNormalAreaTri3D(n, area, p2, p0, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  bool is_lighting = glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_LINES);
  for (size_t itet = 0; itet < aTet.size() / 4; itet++) {
    const unsigned int ino0 = aTet[itet * 4 + 0];
    const unsigned int ino1 = aTet[itet * 4 + 1];
    const unsigned int ino2 = aTet[itet * 4 + 2];
    const unsigned int ino3 = aTet[itet * 4 + 3];
    const double p0[3] = {aXYZ[ino0 * 3 + 0], aXYZ[ino0 * 3 + 1], aXYZ[ino0 * 3 + 2]};
    const double p1[3] = {aXYZ[ino1 * 3 + 0], aXYZ[ino1 * 3 + 1], aXYZ[ino1 * 3 + 2]};
    const double p2[3] = {aXYZ[ino2 * 3 + 0], aXYZ[ino2 * 3 + 1], aXYZ[ino2 * 3 + 2]};
    const double p3[3] = {aXYZ[ino3 * 3 + 0], aXYZ[ino3 * 3 + 1], aXYZ[ino3 * 3 + 2]};
    if (color_glold::IsAbovePlane(p0, org, ncut)) continue;
    if (color_glold::IsAbovePlane(p1, org, ncut)) continue;
    if (color_glold::IsAbovePlane(p2, org, ncut)) continue;
    if (color_glold::IsAbovePlane(p3, org, ncut)) continue;
    // -----
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  //
  /*
   ::glColor3d(0,0,0);
   ::glPointSize(3);
   ::glBegin(GL_POINTS);
   for(unsigned int ino=0;ino<nXYZ_;ino++){
   ::glVertex3dv(pXYZ_+ino*3);
   }
   ::glEnd();
   */
  if (is_lighting) { glEnable(GL_LIGHTING); }
}

// -------------------------

#if defined(_MSC_VER)
#pragma warning( pop )
#endif
