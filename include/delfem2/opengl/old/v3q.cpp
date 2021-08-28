/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/v3q.h"

#include <cassert>

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before opengl headers
#endif

#if defined(__APPLE__) && defined(__MACH__) // mac
#  define GL_SILENCE_DEPRECATION
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

#include "delfem2/vec3.h"
#include "delfem2/quat.h"

#ifndef M_PI
#  define M_PI 3.141592653589793238462643383279502
#endif

//----------------------------------------------------

namespace delfem2 {
namespace opengl {

template<>
DFM2_INLINE void myGlVertex(const CVec3d &v) { ::glVertex3d(v.x, v.y, v.z); }

template<>
DFM2_INLINE void myGlVertex(const CVec3f &v) { ::glVertex3f(v.x, v.y, v.z); }

}
}

// -------------------------

namespace delfem2 {
namespace opengl {

template<>
DFM2_INLINE void myGlTranslate(const CVec3f &v) { ::glTranslatef(v.x, v.y, v.z); }

template<>
DFM2_INLINE void myGlTranslate(const CVec3d &v) { ::glTranslated(v.x, v.y, v.z); }

}
}

// ----------------------------

namespace delfem2 {
namespace opengl {

template<>
DFM2_INLINE void myGlNormal(const CVec3d &n) { ::glNormal3d(n.x, n.y, n.z); }

template<>
DFM2_INLINE void myGlNormal(const CVec3f &n) { ::glNormal3f(n.x, n.y, n.z); }

}
}

// ----------------------------

namespace delfem2 {
namespace opengl {

template<typename REAL>
DFM2_INLINE void myDrawOneTriangle(
    const CVec3<REAL> &v0,
    const CVec3<REAL> &v1,
    const CVec3<REAL> &v2) {
  CVec3<REAL> n;
  UnitNormal(n, v0, v1, v2);
  myGlNormal(n);
  myGlVertex(v0);
  myGlVertex(v1);
  myGlVertex(v2);
}

}
}

// ----------------------------

DFM2_INLINE void delfem2::opengl::myGlNormal(
    const CVec3d &a,
    const CVec3d &b,
    const CVec3d &c) {
  CVec3d n;
  UnitNormal(n, a, b, c);
  ::glNormal3d(n.x, n.y, n.z);
}

DFM2_INLINE void delfem2::opengl::myGlVertex(
    int i,
    const std::vector<CVec3d> &aV) {
  const CVec3d &v = aV[i];
  opengl::myGlVertex(v);
}

DFM2_INLINE void delfem2::opengl::myGlVertex3(
    unsigned int i,
    const std::vector<double> &vec) {
  ::glVertex3d(vec[i * 3 + 0], vec[i * 3 + 1], vec[i * 3 + 2]);
}

DFM2_INLINE void delfem2::opengl::ModelTransformation
    (const CVec3d &dx, const CVec3d &dz, const CVec3d &origin) {
  const CVec3d &dy = Cross(dz, dx);
  const CVec3d &o = origin;
  double A[16];
  A[0] = dx.x;
  A[1] = dx.y;
  A[2] = dx.z;
  A[3] = 0;
  A[4] = dy.x;
  A[5] = dy.y;
  A[6] = dy.z;
  A[7] = 0;
  A[8] = dz.x;
  A[9] = dz.y;
  A[10] = dz.z;
  A[11] = 0;
  A[12] = +o.x;
  A[13] = +o.y;
  A[14] = +o.z;
  A[15] = 1;
  ::glMultMatrixd(A);
}

DFM2_INLINE void delfem2::opengl::ViewTransformation(
    const CVec3d &dx,
    const CVec3d &dz,
    const CVec3d &origin) {
  const CVec3d &dy = Cross(dz, dx);
  CVec3d o(dx.dot(origin), dy.dot(origin), dz.dot(origin));
  double A[16];
  A[0] = dx.x;
  A[1] = dy.x;
  A[2] = dz.x;
  A[3] = 0;
  A[4] = dx.y;
  A[5] = dy.y;
  A[6] = dz.y;
  A[7] = 0;
  A[8] = dx.z;
  A[9] = dy.z;
  A[10] = dz.z;
  A[11] = 0;
  A[12] = -o.x;
  A[13] = -o.y;
  A[14] = -o.z;
  A[15] = 1;
  ::glMultMatrixd(A);
}

//--------------------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawCylinderWire
    (const CVec3d &p0,
     const CVec3d &p1,
     double r) {
  const int ndiv = 16;
  double rdiv = 3.1415 * 2.0 / ndiv;
  CVec3d ez = (p1 - p0).normalized();
  CVec3d ex, ey;
  GetVertical2Vector(ez, ex, ey);
  ::glBegin(GL_LINES);
  for (int idiv = 0; idiv < ndiv; ++idiv) {
    double tA = rdiv * idiv;
    double tB = rdiv * ((idiv + 1) % ndiv);
    CVec3d q0A = p0 + r * ex * sin(tA) + r * ey * cos(tA);
    CVec3d q1A = p1 + r * ex * sin(tA) + r * ey * cos(tA);
    CVec3d q0B = p0 + r * ex * sin(tB) + r * ey * cos(tB);
    CVec3d q1B = p1 + r * ex * sin(tB) + r * ey * cos(tB);
    myGlVertex(q0A);
    myGlVertex(q1A);
    myGlVertex(q0A);
    myGlVertex(q0B);
    myGlVertex(q1A);
    myGlVertex(q1B);
    myGlVertex(p0);
    myGlVertex(q0A);
    myGlVertex(p1);
    myGlVertex(q1A);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawCylinder
    (const CVec3d &p0,
     const CVec3d &p1,
     double r) {
  CVec3d z = (p1 - p0).normalized();
  CVec3d x, y;
  GetVertical2Vector(z, x, y);
  const int ndivt = 32;
  const double dt = 3.1415 * 2.0 / ndivt;
  { // cylinder
    ::glBegin(GL_QUADS);
    for (int idiv = 0; idiv < ndivt; ++idiv) {
      CVec3d n = cos((idiv + 0.5) * dt) * y + sin((idiv + 0.5) * dt) * x;
      myGlNormal(n);
      myGlVertex(p0 + r * sin((idiv + 0) * dt) * x + r * cos((idiv + 0) * dt) * y);
      myGlVertex(p1 + r * sin((idiv + 0) * dt) * x + r * cos((idiv + 0) * dt) * y);
      myGlVertex(p1 + r * sin((idiv + 1) * dt) * x + r * cos((idiv + 1) * dt) * y);
      myGlVertex(p0 + r * sin((idiv + 1) * dt) * x + r * cos((idiv + 1) * dt) * y);
    }
    ::glEnd();
  }
  {
    ::glBegin(GL_TRIANGLES);
    for (int idiv = 0; idiv < ndivt; idiv++) {
      CVec3d v0 = p1 + (r * sin((idiv + 0) * dt)) * x + (r * cos((idiv + 0) * dt)) * y;
      CVec3d v1 = p1 + (r * sin((idiv + 1) * dt)) * x + (r * cos((idiv + 1) * dt)) * y;
      const CVec3d &v2 = p1;
      CVec3d n;
      UnitNormal(n, v1, v0, v2);
      myGlNormal(n);
      myGlVertex(v0);
      myGlVertex(v2);
      myGlVertex(v1);
    }
    ::glEnd();
  }
}

template<typename REAL>
DFM2_INLINE void delfem2::opengl::DrawArrow(
    const CVec3<REAL> &p0,
    const CVec3<REAL> &d,
    int ndivt) {
  using CV3 = CVec3<REAL>;
  CV3 z = d;
  z.normalize();
  CV3 x, y;
  GetVertical2Vector(z, x, y);
  REAL dt = static_cast<REAL>(3.1415 * 2.0) / ndivt;
  REAL r0 = d.norm() * static_cast<REAL>(0.05);
  REAL r1 = d.norm() * static_cast<REAL>(0.10);
  //  double l = d.Length();
  { // cylinder
    ::glBegin(GL_QUADS);
    for (int idiv = 0; idiv < ndivt; idiv++) {
      CV3 n = (REAL) cos((idiv + 0.5) * dt) * y + (REAL) sin((idiv + 0.5) * dt) * x;
      myGlNormal(n);
      CV3 p1 = p0 + d * (REAL) 0.8;
      REAL s0 = (REAL) (r0 * sin((idiv + 0) * dt));
      REAL s1 = (REAL) (r0 * sin((idiv + 1) * dt));
      REAL c0 = (REAL) (r0 * cos((idiv + 0) * dt));
      REAL c1 = (REAL) (r0 * cos((idiv + 1) * dt));
      myGlVertex(p0 + s0 * x + c0 * y);
      myGlVertex(p1 + s0 * x + c0 * y);
      myGlVertex(p1 + s1 * x + c1 * y);
      myGlVertex(p0 + s1 * x + c1 * y);
    }
    ::glEnd();
  }
  { // cone
    ::glBegin(GL_TRIANGLES);
    for (int idiv = 0; idiv < ndivt; idiv++) {
      const CV3 p1 = p0 + d * (REAL) 0.8;
      const REAL s0 = (REAL) (r1 * sin((idiv + 0) * dt));
      const REAL s1 = (REAL) (r1 * sin((idiv + 1) * dt));
      const REAL c0 = (REAL) (r1 * cos((idiv + 0) * dt));
      const REAL c1 = (REAL) (r1 * cos((idiv + 1) * dt));
      const CV3 v0 = p1 + s0 * x + c0 * y;
      const CV3 v1 = p1 + s1 * x + c1 * y;
      const CV3 v2 = p0 + d;
      myDrawOneTriangle(v1, v0, v2);
    }
    ::glEnd();
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::DrawArrow(
    const CVec3d& p0,
    const CVec3d& d,
    int ndivt);
template void delfem2::opengl::DrawArrow(
    const CVec3f& p0,
    const CVec3f& d,
    int ndivt);
#endif

// -----------------------------

template<typename REAL>
DFM2_INLINE void delfem2::opengl::DrawArrowOcta_FaceNrm(
    const CVec3<REAL> &p0,
    const CVec3<REAL> &d,
    REAL rad_ratio,
    REAL node_ratio) {
  using CV3 = CVec3<REAL>;
  CV3 z = d;
  z.normalize();
  CV3 x, y;
  GetVertical2Vector(z, x, y);
  const REAL dt = M_PI * 0.5;
  const REAL r0 = d.norm() * rad_ratio;
  const CV3 p1 = p0 + node_ratio * d;
  const CV3 p2 = p0 + d;
  //
  ::glBegin(GL_TRIANGLES);
  for (int idiv = 0; idiv < 4; idiv++) {
    const REAL s0 = (REAL) (r0 * sin((idiv + 0) * dt));
    const REAL s1 = (REAL) (r0 * sin((idiv + 1) * dt));
    const REAL c0 = (REAL) (r0 * cos((idiv + 0) * dt));
    const REAL c1 = (REAL) (r0 * cos((idiv + 1) * dt));
    const CV3 q0 = p1 + s0 * x + c0 * y;
    const CV3 q1 = p1 + s1 * x + c1 * y;
    myDrawOneTriangle(p0, q0, q1);
    myDrawOneTriangle(p2, q1, q0);
  }
  ::glEnd();
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::DrawArrowOcta_FaceNrm(
    const CVec3d &p0,
    const CVec3d &d,
    double rad_ratio,
    double node_ratio);
template void delfem2::opengl::DrawArrowOcta_FaceNrm(
    const CVec3f &p0,
    const CVec3f &d,
    float rad_ratio,
    float node_ratio);
#endif

template<typename REAL>
DFM2_INLINE void delfem2::opengl::DrawArrowOcta_Edge(
    const delfem2::CVec3<REAL> &p0,
    const delfem2::CVec3<REAL> &d,
    REAL rad_ratio,
    REAL node_ratio) {
  using CV3 = CVec3<REAL>;
  CV3 z = d;
  z.normalize();
  CV3 x, y;
  GetVertical2Vector(z, x, y);
  const REAL dt = M_PI * 0.5;
  const REAL r0 = d.norm() * rad_ratio;
  const CV3 p1 = p0 + node_ratio * d;
  const CV3 p2 = p0 + d;
  //
  ::glBegin(GL_LINES);
  for (int idiv = 0; idiv < 4; idiv++) {
    const REAL s0 = (REAL) (r0 * sin((idiv + 0) * dt));
    const REAL s1 = (REAL) (r0 * sin((idiv + 1) * dt));
    const REAL c0 = (REAL) (r0 * cos((idiv + 0) * dt));
    const REAL c1 = (REAL) (r0 * cos((idiv + 1) * dt));
    const CV3 q0 = p1 + s0 * x + c0 * y;
    const CV3 q1 = p1 + s1 * x + c1 * y;
    myGlVertex(p0);
    myGlVertex(q0);
    myGlVertex(p2);
    myGlVertex(q0);
    myGlVertex(q0);
    myGlVertex(q1);
  }
  ::glEnd();
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::DrawArrowOcta_Edge(
    const delfem2::CVec3d &p0,
    const delfem2::CVec3d &d,
    double rad_ratio,
    double node_ratio);
template void delfem2::opengl::DrawArrowOcta_Edge(
    const delfem2::CVec3f &p0,
    const delfem2::CVec3f &d,
    float rad_ratio,
    float node_ratio);
#endif

// -----------------------------

DFM2_INLINE void delfem2::opengl::DrawCircleArrow(
    const CVec3d &org,
    CVec3d axis,
    double offset) {
  double arrow_width_ratio = 0.1;
  double head_width_ratio = 2.0;
  CVec3d z = -axis;
  z.normalize();
  CVec3d x, y;
  GetVertical2Vector(z, x, y);
  CVec3d p0 = org + offset * z;
  int ndivt = 32;
  double dt = 3.1415 * 2.0 / ndivt;
  double r0 = axis.norm() * arrow_width_ratio;
  double l = axis.norm();
  { // cylinder
    ::glBegin(GL_QUADS);
    for (int idiv = (int) (ndivt * 0.1); idiv < (int) (ndivt * 0.9); idiv++) {
      CVec3d q0 = p0 + l * sin((idiv + 0) * dt) * x + l * cos((idiv + 0) * dt) * y;
      CVec3d q1 = p0 + l * sin((idiv + 1) * dt) * x + l * cos((idiv + 1) * dt) * y;
      CVec3d s0 = sin((idiv + 0) * dt) * x + cos((idiv + 0) * dt) * y;
      CVec3d s1 = sin((idiv + 1) * dt) * x + cos((idiv + 1) * dt) * y;
      for (int jdiv = 0; jdiv < ndivt; jdiv++) {
        CVec3d n = sin((jdiv + 0) * dt) * s0 + cos((jdiv + 0) * dt) * z;
        ::glNormal3d(n.x, n.y, n.z);
        myGlVertex(q0 + r0 * sin((jdiv + 0) * dt) * s0 + r0 * cos((jdiv + 0) * dt) * z);
        myGlVertex(q0 + r0 * sin((jdiv + 1) * dt) * s0 + r0 * cos((jdiv + 1) * dt) * z);
        myGlVertex(q1 + r0 * sin((jdiv + 1) * dt) * s1 + r0 * cos((jdiv + 1) * dt) * z);
        myGlVertex(q1 + r0 * sin((jdiv + 0) * dt) * s1 + r0 * cos((jdiv + 0) * dt) * z);
      }
    }
    ::glEnd();
  }
  { // cone
    double r1 = axis.norm() * head_width_ratio * arrow_width_ratio;
    ::glBegin(GL_TRIANGLES);
    int idiv0 = (int) (ndivt * 0.9 + 1);
    int idiv1 = (int) (ndivt * 1.0);
    CVec3d q0 = p0 + l * sin(idiv0 * dt) * x + l * cos(idiv0 * dt) * y;
    CVec3d q1 = p0 + l * sin(idiv1 * dt) * x + l * cos(idiv1 * dt) * y;
    CVec3d s0 = sin(idiv0 * dt) * x + cos(idiv0 * dt) * y;
    CVec3d s1 = sin(idiv1 * dt) * x + cos(idiv1 * dt) * y;
    for (int jdiv = 0; jdiv < ndivt; jdiv++) {
      CVec3d v0 = q0 + r1 * sin((jdiv + 0) * dt) * s0 + r1 * cos((jdiv + 0) * dt) * z;
      CVec3d v1 = q0 + r1 * sin((jdiv + 1) * dt) * s0 + r1 * cos((jdiv + 1) * dt) * z;
      const CVec3d &v2 = q1;
      CVec3d n;
      UnitNormal(n, v0, v2, v1);
      ::glNormal3d(n.x, n.y, n.z);
      myGlVertex(v0);
      myGlVertex(v2);
      myGlVertex(v1);
    }
    ::glEnd();
  }
}

// --------------------------------------------------------

template<typename REAL>
DFM2_INLINE void delfem2::opengl::DrawCircleWire(
    const CVec3<REAL> &axis,
    const CVec3<REAL> &org,
    REAL r) {
  const REAL pi = static_cast<REAL>(3.1415926535);
  int n = 32;
  REAL dt0 = 2 * pi / n;
  CVec3<REAL> vh, vw;
  GetVertical2Vector(axis, vh, vw);
  ::glBegin(GL_LINE_STRIP);
  for (int i = 0; i < n + 1; i++) {
    CVec3<REAL> p = org + (REAL) (r * sin(dt0 * i)) * vh + (REAL) (r * cos(dt0 * i)) * vw;
    myGlVertex(p);
  }
  ::glEnd();
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::DrawCircleWire(const CVec3f& axis,
                                              const CVec3f& org,
                                              float r);
template void delfem2::opengl::DrawCircleWire(const CVec3d& axis,
                                              const CVec3d& org,
                                              double r);
#endif

// -------------------------

DFM2_INLINE void delfem2::opengl::DrawCircleSolid
    (const CVec3d &axis,
     const CVec3d &org,
     double r) {
  const CVec3d &z = axis;
  CVec3d x, y;
  GetVertical2Vector(z, x, y);
  const int ndivt = 32;
  const double dt = 3.1415 * 2.0 / ndivt;
  {
    ::glBegin(GL_TRIANGLES);
    for (int idiv = 0; idiv < ndivt; idiv++) {
      CVec3d v0 = org + (r * sin((idiv + 0) * dt)) * x + (r * cos((idiv + 0) * dt)) * y;
      CVec3d v1 = org + (r * sin((idiv + 1) * dt)) * x + (r * cos((idiv + 1) * dt)) * y;
      const CVec3d &v2 = org;
      CVec3d n;
      UnitNormal(n, v1, v0, v2);
      myGlNormal(n);
      myGlVertex(v0);
      myGlVertex(v2);
      myGlVertex(v1);
    }
    ::glEnd();
  }
}

DFM2_INLINE void delfem2::opengl::DrawArcSolid
    (const CVec3d &axis,
     const CVec3d &org,
     double ru, // rin
     double rv, // rout
     double rads,
     double rade) {
  const CVec3d &z = axis;
  CVec3d x, y;
  GetVertical2Vector(z, x, y);
  const int ndivt = 32;
  const double dt = (rade - rads) / ndivt;
  {
    ::glBegin(GL_QUADS);
    for (int idiv = 0; idiv < ndivt; idiv++) {
      CVec3d u0 = org + (ru * sin(rads + (idiv + 0) * dt)) * y + (ru * cos(rads + (idiv + 0) * dt)) * x;
      CVec3d u1 = org + (ru * sin(rads + (idiv + 1) * dt)) * y + (ru * cos(rads + (idiv + 1) * dt)) * x;
      CVec3d v0 = org + (rv * sin(rads + (idiv + 0) * dt)) * y + (rv * cos(rads + (idiv + 0) * dt)) * x;
      CVec3d v1 = org + (rv * sin(rads + (idiv + 1) * dt)) * y + (rv * cos(rads + (idiv + 1) * dt)) * x;
      CVec3d n;
      UnitNormal(n, v1, v0, org);
      myGlNormal(n);
      myGlVertex(v0);
      myGlVertex(v1);
      myGlVertex(u1);
      myGlVertex(u0);
    }
    ::glEnd();
  }
}

// ------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawSingleQuad_Edge
    (const CVec3d &p0, const CVec3d &p1, const CVec3d &p2, const CVec3d &p3) {
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINE_LOOP);
  myGlVertex(p0);
  myGlVertex(p1);
  myGlVertex(p2);
  myGlVertex(p3);
  glEnd();
}

// ------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawSingleQuad_FaceNorm
    (const CVec3d &p0, const CVec3d &p1, const CVec3d &p2, const CVec3d &p3) {
  ::glBegin(GL_QUADS);
  {
    CVec3d n0;
    UnitNormal(n0, p0, p1, p3);
    myGlNormal(n0);
    myGlVertex(p0);
  }
  {
    CVec3d n1;
    UnitNormal(n1, p0, p1, p2);
    myGlNormal(n1);
    myGlVertex(p1);
  }
  {
    CVec3d n2;
    UnitNormal(n2, p1, p2, p3);
    myGlNormal(n2);
    myGlVertex(p2);
  }
  {
    CVec3d n3;
    UnitNormal(n3, p2, p3, p0);
    myGlNormal(n3);
    myGlVertex(p3);
  }
  ::glEnd();
}

// ------------------------------------------

DFM2_INLINE void delfem2::opengl::drawPolyLine3D
    (const std::vector<CVec3d> &aP) {
  if (aP.empty()) return;
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip < aP.size() - 1; ip++) {
    unsigned int jp = ip + 1;
    myGlVertex(ip, aP);
    myGlVertex(jp, aP);
  }
  ::glEnd();
  //
  ::glBegin(GL_POINTS);
  for (unsigned int ip = 0; ip < aP.size(); ip++) {
    myGlVertex(ip, aP);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawTriMeshNorm
    (const std::vector<CVec3d> &aP,
     const std::vector<int> &aTri) {
  const int nTri = (int) aTri.size() / 3;
  ::glBegin(GL_TRIANGLES);
  for (int itri = 0; itri < nTri; itri++) {
    const int i0 = aTri[itri * 3 + 0];
    const int i1 = aTri[itri * 3 + 1];
    const int i2 = aTri[itri * 3 + 2];
    const CVec3d &v0 = aP[i0];
    const CVec3d &v1 = aP[i1];
    const CVec3d &v2 = aP[i2];
    const CVec3d &n = Normal(v0, v1, v2).normalized();
    myGlNormal(n);
    myGlVertex(v0);
    myGlVertex(v1);
    myGlVertex(v2);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshTri_Edge
    (const std::vector<CVec3d> &aP,
     const std::vector<unsigned int> &aTri) {
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ////
  /////
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  const int nTri = (int) aTri.size() / 3;
  for (int itri = 0; itri < nTri; ++itri) {
    const int i0 = aTri[itri * 3 + 0];
    const int i1 = aTri[itri * 3 + 1];
    const int i2 = aTri[itri * 3 + 2];
    if (i0 == -1) {
      assert(i1 == -1);
      assert(i2 == -1);
      continue;
    }
    assert(i0 >= 0 && i0 < (int) aP.size());
    assert(i1 >= 0 && i1 < (int) aP.size());
    assert(i2 >= 0 && i2 < (int) aP.size());
    myGlVertex(aP[i0]);
    myGlVertex(aP[i1]);
    myGlVertex(aP[i1]);
    myGlVertex(aP[i2]);
    myGlVertex(aP[i2]);
    myGlVertex(aP[i0]);
  }
  ::glEnd();

  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}

DFM2_INLINE void delfem2::opengl::DrawMeshQuad_Face
    (const std::vector<CVec3d> &aPoint,
     const std::vector<unsigned int> &aQuad) {
  ::glBegin(GL_QUADS);
  for (int iq = 0; iq < (int) aQuad.size() / 4; ++iq) {
    int iv0 = aQuad[iq * 4 + 0];
    int iv1 = aQuad[iq * 4 + 1];
    int iv2 = aQuad[iq * 4 + 2];
    int iv3 = aQuad[iq * 4 + 3];
    {
      CVec3d v01 = aPoint[iv1] - aPoint[iv0];
      CVec3d v12 = aPoint[iv2] - aPoint[iv1];
      CVec3d n = (v01 ^ v12).normalized();
      myGlNormal(n);
    }
    myGlVertex(aPoint[iv0]);
    myGlVertex(aPoint[iv1]);
    myGlVertex(aPoint[iv2]);
    myGlVertex(aPoint[iv3]);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawPoint3D
    (const std::vector<CVec3d> &aPoint) {
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_POINTS);
  for (const auto &p : aPoint) {
    myGlVertex(p);
  }
  ::glEnd();
}

// ------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawQuad3D_Edge
    (const std::vector<CVec3d> &aPoint,
     const std::vector<unsigned int> &aQuad) {
  ::glBegin(GL_LINES);
  for (int iq = 0; iq < (int) aQuad.size() / 4; ++iq) {
    int iv0 = aQuad[iq * 4 + 0];
    int iv1 = aQuad[iq * 4 + 1];
    int iv2 = aQuad[iq * 4 + 2];
    int iv3 = aQuad[iq * 4 + 3];
    myGlVertex(aPoint[iv0]);
    myGlVertex(aPoint[iv1]);
    myGlVertex(aPoint[iv1]);
    myGlVertex(aPoint[iv2]);
    myGlVertex(aPoint[iv2]);
    myGlVertex(aPoint[iv3]);
    myGlVertex(aPoint[iv3]);
    myGlVertex(aPoint[iv0]);
  }
  ::glEnd();
}

// ------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawSingleHex_Edge
    (const CVec3d &p0, const CVec3d &p1, const CVec3d &p2, const CVec3d &p3,
     const CVec3d &p4, const CVec3d &p5, const CVec3d &p6, const CVec3d &p7) {
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINE_LOOP);
  myGlVertex(p0);
  myGlVertex(p1);
  myGlVertex(p2);
  myGlVertex(p3);
  glEnd();
  ::glBegin(GL_LINE_LOOP);
  myGlVertex(p4);
  myGlVertex(p5);
  myGlVertex(p6);
  myGlVertex(p7);
  glEnd();
  ::glBegin(GL_LINES);
  myGlVertex(p0);
  myGlVertex(p4);
  myGlVertex(p1);
  myGlVertex(p5);
  myGlVertex(p2);
  myGlVertex(p6);
  myGlVertex(p3);
  myGlVertex(p7);
  ::glEnd();
}

// ------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawGrid2D
    (int ndivx, int ndivy,
     const CVec3d &ex, const CVec3d &ey, const CVec3d &org) {
  const CVec3d &p00 = org;
  const CVec3d &p10 = org + ex * (double) ndivx;
  const CVec3d &p01 = org + ey * (double) ndivy;
  ::glBegin(GL_LINES);
  for (int ix = 0; ix < ndivx + 1; ++ix) {
    const CVec3d &dx = (double) ix * ex;
    myGlVertex(p00 + dx);
    myGlVertex(p01 + dx);
  }
  for (int iy = 0; iy < ndivy + 1; ++iy) {
    const CVec3d &dy = (double) iy * ey;
    myGlVertex(p00 + dy);
    myGlVertex(p10 + dy);
  }
  ::glEnd();
}

// ------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawGridOutside
    (int ndivx, int ndivy, int ndivz,
     double elen,
     const CVec3d &org) {
  DrawGrid2D(ndivx, ndivy, CVec3d(elen, 0, 0), CVec3d(0, elen, 0), org);
  DrawGrid2D(ndivx, ndivy, CVec3d(elen, 0, 0), CVec3d(0, elen, 0), org + CVec3d(0, 0, elen * ndivz));
  DrawGrid2D(ndivy, ndivz, CVec3d(0, elen, 0), CVec3d(0, 0, elen), org);
  DrawGrid2D(ndivy, ndivz, CVec3d(0, elen, 0), CVec3d(0, 0, elen), org + CVec3d(elen * ndivx, 0, 0));
  DrawGrid2D(ndivz, ndivx, CVec3d(0, 0, elen), CVec3d(elen, 0, 0), org);
  DrawGrid2D(ndivz, ndivx, CVec3d(0, 0, elen), CVec3d(elen, 0, 0), org + CVec3d(0, elen * ndivy, 0));
}

// -----------------------------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawHandlerRotation_Mat4
    (const double Mat[16],
     double size,
     int ielem_picked) {
  ::glDisable(GL_LIGHTING);
  {
    if (ielem_picked == 0) { ::glColor3d(1, 1, 0); } else { ::glColor3d(1, 0, 0); }
    const CVec3d &ax = Mat4Vec(Mat, CVec3d(1, 0, 0)).normalized();
    const CVec3d pos(Mat[3], Mat[7], Mat[11]);
    opengl::DrawCircleWire(ax, pos, size);
  }
  {
    if (ielem_picked == 1) { ::glColor3d(1, 1, 0); } else { ::glColor3d(0, 1, 0); }
    const CVec3d &ay = Mat4Vec(Mat, CVec3d(0, 1, 0)).normalized();
    const CVec3d pos(Mat[3], Mat[7], Mat[11]);
    opengl::DrawCircleWire(ay, pos, size);
  }
  {
    if (ielem_picked == 2) { ::glColor3d(1, 1, 0); } else { ::glColor3d(0, 0, 1); }
    const CVec3d &az = Mat4Vec(Mat, CVec3d(0, 0, 1)).normalized();
    const CVec3d pos(Mat[3], Mat[7], Mat[11]);
    opengl::DrawCircleWire(az, pos, size);
  }
}

// ----------------------------------

namespace delfem2 {
namespace opengl {

template<>
DFM2_INLINE void MyGlMultMat(const CMat4<double> &m) {
  glMultMatrixd(m.mat);
}

template<>
DFM2_INLINE void MyGlMultMat(const CMat4<float> &m) {
  glMultMatrixf(m.mat);
}

}
}

// above: gizmo
// -------------------------------------------------------
// below: quaternion coord

DFM2_INLINE void delfem2::opengl::Draw_QuaternionsCoordinateAxes
    (const std::vector<double> &aXYZ1,
     const std::vector<double> &aQuat,
     double l) {
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(2);
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip < aXYZ1.size() / 3; ++ip) {
    const double *p = aXYZ1.data() + ip * 3;
    {
      double ex0[3] = {1, 0, 0};
      double ex[3];
      QuatVec(ex, aQuat.data() + ip * 4, ex0);
      ::glColor3d(1, 0, 0);
      ::glVertex3dv(p);
      ::glVertex3d(p[0] + l * ex[0], p[1] + l * ex[1], p[2] + l * ex[2]);
    }
    {
      double ey0[3] = {0, 1, 0};
      double ey[3];
      QuatVec(ey, aQuat.data() + ip * 4, ey0);
      ::glColor3d(0, 1, 0);
      ::glVertex3dv(p);
      ::glVertex3d(p[0] + l * ey[0], p[1] + l * ey[1], p[2] + l * ey[2]);
    }
    {
      double ez0[3] = {0, 0, 1};
      double ez[3];
      QuatVec(ez, aQuat.data() + ip * 4, ez0);
      ::glColor3d(0, 0, 1);
      ::glVertex3dv(p);
      ::glVertex3d(p[0] + l * ez[0], p[1] + l * ez[1], p[2] + l * ez[2]);
    }
  }
  ::glEnd();
}

