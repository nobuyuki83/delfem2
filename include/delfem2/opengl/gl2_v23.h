/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GL2_V23_H
#define DFM2_GL2_V23_H

#include <vector>
#include "delfem2/vec2.h"
#include "delfem2/vec3.h"

namespace delfem2{
namespace opengl
{

void myGlVertex(const CVector3& v);
void myGlTranslate(const CVector3& v);
void myGlNormal(const CVector3& n);
void myGlNormal(const CVector3& a, const CVector3& b, const CVector3& c);
void myGlVertex(int i, const std::vector<CVector3>& aV);
void myGlVertex(int i, const std::vector<double>& vec);
void ModelTransformation(const CVector3& dx, const CVector3& dz, const CVector3& origin);
void ViewTransformation(const CVector3& dx, const CVector3& dz, const CVector3& origin);

// -------------------------------------

void DrawArcSolid(const CVector3& axis,
                  const CVector3& org,
                  double ru, // rin
                  double rv, // rout
                  double rads,
                  double rade);
void DrawArrow(const CVector3& p0,
               const CVector3& d,
               int ndivt=16);
void DrawCircleArrow(const CVector3& org, CVector3 axis, double offset);
void DrawCylinder(const CVector3& p0,
                  const CVector3& p1,
                  double r);
void DrawCylinderWire(const CVector3& p0,
                      const CVector3& p1,
                      double r);
void DrawSingleQuad_FaceNorm(const CVector3& p0,
                             const CVector3& p1,
                             const CVector3& p2,
                             const CVector3& p3);
void DrawSingleQuad_Edge(const CVector3& p0,
                         const CVector3& p1,
                         const CVector3& p2,
                         const CVector3& p3);
void DrawSingleHex_Edge(const CVector3& p0,
                        const CVector3& p1,
                        const CVector3& p2,
                        const CVector3& p3,
                        const CVector3& p4,
                        const CVector3& p5,
                        const CVector3& p6,
                        const CVector3& p7);
void drawPolyLine3D(const std::vector<CVector3>& aP);

void DrawCircleWire(const CVector3& axis,
                    const CVector3& org,
                    double r);
void DrawCircleSolid(const CVector3& axis,
                     const CVector3& org,
                     double r);
void DrawGrid2D(int ndivx, int ndivy,
                const CVector3& ex, const CVector3& ey, const CVector3& org);
void DrawGridOutside(int ndivx, int ndivy, int ndivz,
                     double elen,
                     const CVector3& org);
  

// ------------
// mesh from here
void DrawPoint3D(const std::vector<CVector3>& aPoint);
void drawPolyLine(const std::vector<CVector2>& aP);
void DrawMeshQuad_Face(const std::vector<CVector3>& aPoint,
                       const std::vector<unsigned int>& aQuad);
void DrawMeshTri_Edge(const std::vector<CVector3>& aP,
                      const std::vector<unsigned int>& aTri);
void DrawTriMeshNorm(const std::vector<CVector3>& aP,
                     const std::vector<int>& aTri);
void DrawMeshTri_Edge(const std::vector<CVector3>& aP,
                      const std::vector<unsigned int>& aTri);
void DrawQuad3D_Edge(const std::vector<CVector3>& aPoint,
                     const std::vector<unsigned int>& aQuad);

// -----------
// Handler

void DrawAxisHandler(double s, const CVector3& p);
void DrawHandlerRotation_PosQuat(const CVector3& pos, const double quat[4],
                                 double size, int ielem_picked);
void DrawHandlerRotation_Mat4(const double Mat[16],
                              double size, int ielem_picked);

// ----------------
// vec2 starts here

void myGlVertex(int i, const std::vector<CVector2>& aP);

void myGlVertex(const CVector2& v);

void drawPolyLine2D(const std::vector<CVector2>& aP);

void Draw_MeshTri(const std::vector<CVector2>& aP,
                  const std::vector<unsigned int>& aTri);

void Draw_MeshTri_Edge(const std::vector<CVector2>& aP,
                       const std::vector<unsigned int>& aTri);
  
} // end namespace opengl
} // end namespace delfem2

#endif
