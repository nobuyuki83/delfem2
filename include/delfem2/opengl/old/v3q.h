/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_OLD_V3Q_H
#define DFM2_OPENGL_OLD_V3Q_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"

namespace delfem2 {
namespace opengl {

// ------------------------------------------------------------------------------------
// vec3 starts here

template<typename REAL>
DFM2_INLINE void myGlVertex(const CVec3<REAL> &v);

template<typename REAL>
DFM2_INLINE void myGlTranslate(const CVec3<REAL> &v);

template<typename REAL>
DFM2_INLINE void myGlNormal(const CVec3<REAL> &n);

DFM2_INLINE void myGlNormal(const CVec3d &a, const CVec3d &b, const CVec3d &c);
DFM2_INLINE void myGlVertex(int i, const std::vector<CVec3d> &aV);
DFM2_INLINE void myGlVertex3(unsigned int i, const std::vector<double> &vec);
DFM2_INLINE void ModelTransformation(const CVec3d &dx, const CVec3d &dz, const CVec3d &origin);
DFM2_INLINE void ViewTransformation(const CVec3d &dx, const CVec3d &dz, const CVec3d &origin);

// -------------------------------------

DFM2_INLINE void DrawArcSolid(
    const CVec3d &axis,
    const CVec3d &org,
    double ru, // rin
    double rv, // rout
    double rads,
    double rade);

template<typename REAL>
DFM2_INLINE void DrawArrow(
    const CVec3<REAL> &p0,
    const CVec3<REAL> &d,
    int ndivt = 16);

template<typename REAL>
DFM2_INLINE void DrawArrowOcta_FaceNrm(
    const CVec3<REAL> &p,
    const CVec3<REAL> &d,
    REAL rad_ratio,
    REAL node_ratio);

template<typename REAL>
DFM2_INLINE void DrawArrowOcta_Edge(
    const delfem2::CVec3<REAL> &p,
    const delfem2::CVec3<REAL> &d,
    REAL rad_ratio,
    REAL node_ratio);

DFM2_INLINE void DrawCircleArrow(
    const CVec3d &org,
    CVec3d axis,
    double offset);

DFM2_INLINE void DrawCylinder(
    const CVec3d &p0,
    const CVec3d &p1,
    double r);

DFM2_INLINE void DrawCylinderWire(
    const CVec3d &p0,
    const CVec3d &p1,
    double r);

DFM2_INLINE void DrawSingleQuad_FaceNorm(
    const CVec3d &p0,
    const CVec3d &p1,
    const CVec3d &p2,
    const CVec3d &p3);

DFM2_INLINE void DrawSingleQuad_Edge(
    const CVec3d &p0,
    const CVec3d &p1,
    const CVec3d &p2,
    const CVec3d &p3);

DFM2_INLINE void DrawSingleHex_Edge(
    const CVec3d &p0,
    const CVec3d &p1,
    const CVec3d &p2,
    const CVec3d &p3,
    const CVec3d &p4,
    const CVec3d &p5,
    const CVec3d &p6,
    const CVec3d &p7);

DFM2_INLINE void drawPolyLine3D(
    const std::vector<CVec3d> &aP);

template<typename REAL>
DFM2_INLINE void DrawCircleWire(
    const CVec3<REAL> &axis,
    const CVec3<REAL> &org,
    REAL r);

DFM2_INLINE void DrawCircleSolid(
    const CVec3d &axis,
    const CVec3d &org,
    double r);

DFM2_INLINE void DrawGrid2D(
    int ndivx,
    int ndivy,
    const CVec3d &ex,
    const CVec3d &ey,
    const CVec3d &org);

void DrawGridOutside(
    int ndivx,
    int ndivy,
    int ndivz,
    double elen,
    const CVec3d &org);

// ------------
// mesh from here

DFM2_INLINE void DrawPoint3D(
    const std::vector<CVec3d> &aPoint);

DFM2_INLINE void DrawMeshQuad_Face(
    const std::vector<CVec3d> &aPoint,
    const std::vector<unsigned int> &aQuad);

DFM2_INLINE void DrawMeshTri_Edge(
    const std::vector<CVec3d> &aP,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void DrawTriMeshNorm(
    const std::vector<CVec3d> &aP,
    const std::vector<int> &aTri);

DFM2_INLINE void DrawMeshTri_Edge(
    const std::vector<CVec3d> &aP,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void DrawQuad3D_Edge(
    const std::vector<CVec3d> &aPoint,
    const std::vector<unsigned int> &aQuad);

// CVec3
// --------------------------------------------------------------------
// CMat4

template<typename REAL>
DFM2_INLINE void MyGlMultMat(
    const delfem2::CMat4<REAL> &m);

DFM2_INLINE void DrawHandlerRotation_Mat4(
    const double Mat[16],
    double size, int ielem_picked);

// CMat4
// --------------------------------------------------------------------
// CQuaternion

DFM2_INLINE void Draw_QuaternionsCoordinateAxes(
    const std::vector<double> &aXYZ1,
    const std::vector<double> &aQuat,
    double l);

} // end namespace opengl
} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/old/v3q.cpp"
#endif

#endif
