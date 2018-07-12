#ifndef UTILITY_GLVECTOR_H
#define UTILITY_GLVECTOR_H

#include <vector>

#if defined(__APPLE__) && defined(__MACH__) // mac
  #include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
  #include <GL/glu.h>
#elif defined(WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else // linux
  #include <GL/gl.h>
#endif

// eventually, I would like to split vec2 and vec3 dependent part.
#include "vec2.h"
#include "vec3.h"
#include "quat.h"

void myGlVertex(const CVector3& v);
void myGlTranslate(const CVector3& v);
void myGlNormal(const CVector3& n);
void myGlNormal(const CVector3& a, const CVector3& b, const CVector3& c);
void myGlVertex(int i, const std::vector<CVector3>& aV);
void myGlVertex(int i, const std::vector<double>& vec);

///////////////////////////////////////////////////

void DrawArcSolid(const CVector3& axis,
                  const CVector3& org,
                  double ru, // rin
                  double rv, // rout
                  double rads,
                  double rade);
void DrawArrow(const CVector3& p0,
               const CVector3& d,
               int ndivt=16);
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




////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// vec2 starts here

inline void myGlVertex(int i, const std::vector<CVector2>& aP)
{
  ::glVertex3d(aP[i].x, aP[i].y, +0.0);
}

inline void myGlVertex(const CVector2& v)
{
  ::glVertex2d(v.x,v.y);
}

CVector2 screenXYProjection(const CVector3& v,
                            const float* mMV, const float* mPj);

void drawPolyLine2D(const std::vector<CVector2>& aP);



//////////////////////////////////////////////////
// Quaterion from here

CVector3 QuatVec(const double quat[4], const CVector3& v0);
CVector3 QuatConjVec(const double quat[4], const CVector3& v0);

/////////////////
// Handler

void DrawAxisHandler(double s, const CVector3& p);
void DrawHandlerRotation(const CVector3& pos, const double quat[4], double size, int ielem_picked);
int PickHandlerRotation(const CVector3& src, const CVector3& dir,
                        const CVector3& pos, const double quat[4], double rad,
                        double tol);
bool DragHandlerRot(double quat[4], int ielem,
                    const CVector2& sp0, const CVector2& sp1,
                    const CVector3& pos,
                    const float mMV[16], const float mPj[16]);
CVector3 drag_AxisHandler(const CVector2& sp0,
                          const CVector2& sp1,
                          const CVector3& p,
                          const CVector3& axis,
                          double len,
                          const float* mMV,
                          const float* mPj);
bool isPickPoint(const CVector2& sp,
                 const CVector3& p,
                 const float* mMV,
                 const float* mPj,
                 double pick_tol);
bool isPick_AxisHandler(const CVector2& sp,
                        const CVector3& p,
                        const CVector3& axis,
                        double len,
                        const float* mMV,
                        const float* mPj,
                        double pick_tol);
bool isPickQuad(const CVector3& p0,
                const CVector3& p1,
                const CVector3& p2,
                const CVector3& p3,
                const CVector2& sp,
                const CVector3& pick_dir,
                const float mMV[16],
                const float mPj[16],
                double eps);

bool isPickCircle(const CVector3& axis,
                  const CVector3& org,
                  double rad,
                  const CVector3& src,
                  const CVector3& dir,
                  double pick_tol);
double DragCircle(const CVector2& sp0,
                  const CVector2& sp1,
                  const CVector3& p,
                  const CVector3& axis,
                  const float* mMV,
                  const float* mPj);


#endif
