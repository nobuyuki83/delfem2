/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cstdio>
#include <cassert>

#if defined(__APPLE__) && defined(__MACH__) // mac
#include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#else // linux
#include <GL/gl.h>
#endif

#include "delfem2/quat.h"
#include "delfem2/gl_v23q.h"

void myGlVertex(const CVector3& v)
{
  ::glVertex3d(v.x,v.y,v.z);
}

void myGlTranslate(const CVector3& v)
{
  ::glTranslated(v.x,v.y,v.z);
}

void myGlNormal(const CVector3& n)
{
  ::glNormal3d(n.x,n.y,n.z);
}

void myGlNormal(const CVector3& a, const CVector3& b, const CVector3& c)
{
  CVector3 n; UnitNormal(n, a, b, c);
  ::glNormal3d(n.x,n.y,n.z);
}

void myGlVertex(int i, const std::vector<CVector3>& aV)
{
  const CVector3& v = aV[i];
  ::myGlVertex(v);
}

void myGlVertex(int i, const std::vector<double>& vec)
{
  ::glVertex3d(vec[i*2], vec[i*2+1], +0.0);
}

void ModelTransformation(const CVector3& dx, const CVector3& dz, const CVector3& origin)
{
  const CVector3& dy = Cross(dz,dx);
  const CVector3& o = origin;
  double A[16];
  A[ 0] = dx.x;  A[ 1] = dx.y;  A[ 2] = dx.z;  A[ 3] = 0;
  A[ 4] = dy.x;  A[ 5] = dy.y;  A[ 6] = dy.z;  A[ 7] = 0;
  A[ 8] = dz.x;  A[ 9] = dz.y;  A[10] = dz.z;  A[11] = 0;
  A[12] = +o.x;  A[13] =  +o.y; A[14] = +o.z;  A[15] = 1;
  ::glMultMatrixd(A);
}

void ViewTransformation(const CVector3& dx, const CVector3& dz, const CVector3& origin)
{
  const CVector3& dy = Cross(dz,dx);
  CVector3 o(dx*origin,dy*origin,dz*origin);
  double A[16];
  A[ 0] = dx.x;  A[ 1] = dy.x;  A[ 2] = dz.x;  A[ 3] = 0;
  A[ 4] = dx.y;  A[ 5] = dy.y;  A[ 6] = dz.y;  A[ 7] = 0;
  A[ 8] = dx.z;  A[ 9] = dy.z;  A[10] = dz.z;  A[11] = 0;
  A[12] = -o.x;  A[13] = -o.y;  A[14] = -o.z;  A[15] = 1;
  ::glMultMatrixd(A);
}

///////////////////////////////////////////////

CVector2 screenXYProjection
(const CVector3& v,
 const float* mMV,
 const float* mPj)
{
  CVector3 sp0 = screenProjection(v,mMV,mPj);
  return CVector2(sp0.x,sp0.y);
}

/////////////

void DrawCylinderWire
(const CVector3& p0,
 const CVector3& p1,
 double r)
{
  const int ndiv = 16;
  double rdiv = 3.1415*2.0/ndiv;
  CVector3 ez = (p1-p0).Normalize();
  CVector3 ex,ey; GetVertical2Vector(ez, ex, ey);
  ::glBegin(GL_LINES);
  for(int idiv=0;idiv<ndiv;++idiv){
    double tA = rdiv*idiv;
    double tB = rdiv*((idiv+1)%ndiv);
    CVector3 q0A = p0+r*ex*sin(tA)+r*ey*cos(tA);
    CVector3 q1A = p1+r*ex*sin(tA)+r*ey*cos(tA);
    CVector3 q0B = p0+r*ex*sin(tB)+r*ey*cos(tB);
    CVector3 q1B = p1+r*ex*sin(tB)+r*ey*cos(tB);
    myGlVertex(q0A); myGlVertex(q1A);
    myGlVertex(q0A); myGlVertex(q0B);
    myGlVertex(q1A); myGlVertex(q1B);
    myGlVertex(p0); myGlVertex(q0A);
    myGlVertex(p1); myGlVertex(q1A);
  }
  ::glEnd();
}

void DrawCylinder
(const CVector3& p0,
 const CVector3& p1,
 double r)
{
  CVector3 z = (p1-p0).Normalize();
  CVector3 x, y; GetVertical2Vector(z, x, y);
  const int ndivt = 32;
  const double dt = 3.1415*2.0/ndivt;
  { // cylinder
    ::glBegin(GL_QUADS);
    for (int idiv=0; idiv<ndivt; ++idiv){
      CVector3 n = cos((idiv+0.5)*dt)*y+sin((idiv+0.5)*dt)*x;
      ::myGlNormal(n);
      ::myGlVertex(p0+r*sin((idiv+0)*dt)*x+r*cos((idiv+0)*dt)*y);
      ::myGlVertex(p1+r*sin((idiv+0)*dt)*x+r*cos((idiv+0)*dt)*y);
      ::myGlVertex(p1+r*sin((idiv+1)*dt)*x+r*cos((idiv+1)*dt)*y);
      ::myGlVertex(p0+r*sin((idiv+1)*dt)*x+r*cos((idiv+1)*dt)*y);
    }
    ::glEnd();
  }
  {
    ::glBegin(GL_TRIANGLES);
    for (int idiv = 0; idiv<ndivt; idiv++){
      CVector3 v0 = p1+(r*sin((idiv+0)*dt))*x+(r*cos((idiv+0)*dt))*y;
      CVector3 v1 = p1+(r*sin((idiv+1)*dt))*x+(r*cos((idiv+1)*dt))*y;
      CVector3 v2 = p1;
      CVector3 n; UnitNormal(n, v1, v0, v2);
      ::myGlNormal(n);
      ::myGlVertex(v0);
      ::myGlVertex(v2);
      ::myGlVertex(v1);
    }
    ::glEnd();
  }
}


void DrawArrow
(const CVector3& p0,
 const CVector3& d,
 int ndivt)
{
  CVector3 z = d; z.SetNormalizedVector();
  CVector3 x,y; GetVertical2Vector(z,x,y);
  double dt = 3.1415*2.0 / ndivt;
  double r0 = d.Length()*0.05;
  double r1 = d.Length()*0.10;
  //  double l = d.Length();
  { // cylinder
    ::glBegin(GL_QUADS);
    for(int idiv=0;idiv<ndivt;idiv++){
      CVector3 n = cos((idiv+0.5)*dt)*y + sin((idiv+0.5)*dt)*x;
      ::myGlNormal(n);
      ::myGlVertex(p0      +r0*sin((idiv+0)*dt)*x+r0*cos((idiv+0)*dt)*y);
      ::myGlVertex(p0+d*0.8+r0*sin((idiv+0)*dt)*x+r0*cos((idiv+0)*dt)*y);
      ::myGlVertex(p0+d*0.8+r0*sin((idiv+1)*dt)*x+r0*cos((idiv+1)*dt)*y);
      ::myGlVertex(p0      +r0*sin((idiv+1)*dt)*x+r0*cos((idiv+1)*dt)*y);
    }
    ::glEnd();
  }
  { // cone
    ::glBegin(GL_TRIANGLES);
    for(int idiv=0;idiv<ndivt;idiv++){
      CVector3 v0 = p0+d*0.8 + (r1*sin((idiv+0)*dt))*x + (r1*cos((idiv+0)*dt))*y;
      CVector3 v1 = p0+d*0.8 + (r1*sin((idiv+1)*dt))*x + (r1*cos((idiv+1)*dt))*y;
      CVector3 v2 = p0+d;
      CVector3 n; UnitNormal(n, v1, v0, v2);
      ::myGlNormal(n);
      ::myGlVertex(v0);
      ::myGlVertex(v2);
      ::myGlVertex(v1);
    }
    ::glEnd();
  }
}

void DrawCircleArrow
(CVector3 org, CVector3 axis, double offset)
{
  double arrow_width_ratio = 0.1;
  double head_width_ratio = 2.0;
  CVector3 z = -axis; z.SetNormalizedVector();
  CVector3 x,y; GetVertical2Vector(z,x,y);
  CVector3 p0 = org+offset*z;
  int ndivt = 32;
  double dt = 3.1415*2.0 / ndivt;
  double r0 = axis.Length()*arrow_width_ratio;
  double l = axis.Length();
  { // cylinder
    ::glBegin(GL_QUADS);
    for(int idiv=(int)(ndivt*0.1);idiv<(int)(ndivt*0.9);idiv++){
      CVector3 q0 = p0 + l*sin((idiv+0)*dt)*x + l*cos((idiv+0)*dt)*y;
      CVector3 q1 = p0 + l*sin((idiv+1)*dt)*x + l*cos((idiv+1)*dt)*y;
      CVector3 s0 = sin((idiv+0)*dt)*x + cos((idiv+0)*dt)*y;
      CVector3 s1 = sin((idiv+1)*dt)*x + cos((idiv+1)*dt)*y;
      for(int jdiv=0;jdiv<ndivt;jdiv++){
        CVector3 n = sin((jdiv+0)*dt)*s0 + cos((jdiv+0)*dt)*z;
        ::glNormal3d(n.x,n.y,n.z);
        ::myGlVertex(q0 + r0*sin((jdiv+0)*dt)*s0 + r0*cos((jdiv+0)*dt)*z);
        ::myGlVertex(q0 + r0*sin((jdiv+1)*dt)*s0 + r0*cos((jdiv+1)*dt)*z);
        ::myGlVertex(q1 + r0*sin((jdiv+1)*dt)*s1 + r0*cos((jdiv+1)*dt)*z);
        ::myGlVertex(q1 + r0*sin((jdiv+0)*dt)*s1 + r0*cos((jdiv+0)*dt)*z);
      }
    }
    ::glEnd();
  }
  { // cone
    double r1 = axis.Length()*head_width_ratio*arrow_width_ratio;
    ::glBegin(GL_TRIANGLES);
    int idiv0 = (int)(ndivt*0.9+1);
    int idiv1 = (int)(ndivt*1.0);
    CVector3 q0 = p0 + l*sin(idiv0*dt)*x + l*cos(idiv0*dt)*y;
    CVector3 q1 = p0 + l*sin(idiv1*dt)*x + l*cos(idiv1*dt)*y;
    CVector3 s0 = sin(idiv0*dt)*x + cos(idiv0*dt)*y;
    CVector3 s1 = sin(idiv1*dt)*x + cos(idiv1*dt)*y;
    for(int jdiv=0;jdiv<ndivt;jdiv++){
      CVector3 v0 = q0 + r1*sin((jdiv+0)*dt)*s0 + r1*cos((jdiv+0)*dt)*z;
      CVector3 v1 = q0 + r1*sin((jdiv+1)*dt)*s0 + r1*cos((jdiv+1)*dt)*z;
      CVector3 v2 = q1;
      CVector3 n; UnitNormal(n, v0, v2, v1);
      ::glNormal3d(n.x,n.y,n.z);
      ::myGlVertex(v0);
      ::myGlVertex(v2);
      ::myGlVertex(v1);
    }
    ::glEnd();
  }
}



//////////////////////////////////////////////////////

void DrawCircleWire
(const CVector3& axis,
 const CVector3& org,
 double r)
{
  const double pi = 3.1415926535;
  int n = 32; double dt = 2*pi/n;
  CVector3 h,v; GetVertical2Vector(axis, h, v);
  ::glBegin(GL_LINE_STRIP);
  for(int i=0;i<n+1;i++) {
    CVector3 p  = org + (r*sin(dt*i))*h + (r*cos(dt*i))*v;
    myGlVertex(p);
  }
  ::glEnd();
}

void DrawCircleSolid
(const CVector3& axis,
 const CVector3& org,
 double r)
{
  CVector3 z = axis;
  CVector3 x, y; GetVertical2Vector(z, x, y);
  const int ndivt = 32;
  const double dt = 3.1415*2.0/ndivt;
  {
    ::glBegin(GL_TRIANGLES);
    for (int idiv = 0; idiv<ndivt; idiv++){
      CVector3 v0 = org+(r*sin((idiv+0)*dt))*x+(r*cos((idiv+0)*dt))*y;
      CVector3 v1 = org+(r*sin((idiv+1)*dt))*x+(r*cos((idiv+1)*dt))*y;
      CVector3 v2 = org;
      CVector3 n; UnitNormal(n, v1, v0, v2);
      ::myGlNormal(n);
      ::myGlVertex(v0);
      ::myGlVertex(v2);
      ::myGlVertex(v1);
    }
    ::glEnd();
  }
}

void DrawArcSolid
(const CVector3& axis,
 const CVector3& org,
 double ru, // rin
 double rv, // rout
 double rads,
 double rade)
{
  CVector3 z = axis;
  CVector3 x, y; GetVertical2Vector(z, x, y);
  const int ndivt = 32;
  const double dt = (rade-rads)/ndivt;
  {
    ::glBegin(GL_QUADS);
    for (int idiv = 0; idiv<ndivt; idiv++){
      CVector3 u0 = org+(ru*sin(rads+(idiv+0)*dt))*y+(ru*cos(rads+(idiv+0)*dt))*x;
      CVector3 u1 = org+(ru*sin(rads+(idiv+1)*dt))*y+(ru*cos(rads+(idiv+1)*dt))*x;
      CVector3 v0 = org+(rv*sin(rads+(idiv+0)*dt))*y+(rv*cos(rads+(idiv+0)*dt))*x;
      CVector3 v1 = org+(rv*sin(rads+(idiv+1)*dt))*y+(rv*cos(rads+(idiv+1)*dt))*x;
      CVector3 n; UnitNormal(n, v1, v0, org);
      ::myGlNormal(n);
      ::myGlVertex(v0);
      ::myGlVertex(v1);
      ::myGlVertex(u1);
      ::myGlVertex(u0);
    }
    ::glEnd();
  }
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void DrawSingleQuad_Edge
(const CVector3& p0, const CVector3& p1, const CVector3& p2, const CVector3& p3)
{
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINE_LOOP);
  myGlVertex(p0);
  myGlVertex(p1);
  myGlVertex(p2);
  myGlVertex(p3);
  glEnd();
  
}

void DrawSingleQuad_FaceNorm
(const CVector3& p0, const CVector3& p1, const CVector3& p2, const CVector3& p3)
{
  ::glBegin(GL_QUADS);
  {
    CVector3 n0; UnitNormal(n0,  p0, p1, p3);
    myGlNormal(n0);
    myGlVertex(p0);
  }
  {
    CVector3 n1; UnitNormal(n1,  p0, p1, p2);
    myGlNormal(n1);
    myGlVertex(p1);
  }
  {
    CVector3 n2; UnitNormal(n2,  p1, p2, p3);
    myGlNormal(n2);
    myGlVertex(p2);
  }
  {
    CVector3 n3; UnitNormal(n3,  p2, p3, p0);
    myGlNormal(n3);
    myGlVertex(p3);
  }
  ::glEnd();
}


void drawPolyLine
(const std::vector<CVector2>& aP)
{
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip<aP.size()-1; ip++){
    unsigned int jp = ip+1;
    myGlVertex(ip,aP);
    myGlVertex(jp,aP);
  }
  ::glEnd();
  
  ////
  ::glBegin(GL_POINTS);
  for (unsigned int ip = 0; ip<aP.size(); ip++){
    myGlVertex(ip,aP);
  }
  ::glEnd();
}

void drawPolyLine3D(const std::vector<CVector3>& aP)
{
  if( aP.empty() ) return;
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip<aP.size()-1; ip++){
    unsigned int jp = ip+1;
    myGlVertex(ip,aP);
    myGlVertex(jp,aP);
  }
  ::glEnd();
  
  ////
  ::glBegin(GL_POINTS);
  for (unsigned int ip = 0; ip<aP.size(); ip++){
    myGlVertex(ip,aP);
  }
  ::glEnd();
}

void drawPolyLine2D(const std::vector<CVector2>& aP)
{
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip<aP.size()-1; ip++){
    unsigned int jp = ip+1;
    myGlVertex(ip,aP);
    myGlVertex(jp,aP);
  }
  ::glEnd();
  
  ////
  ::glBegin(GL_POINTS);
  for (unsigned int ip = 0; ip<aP.size(); ip++){
    myGlVertex(ip,aP);
  }
  ::glEnd();
}

void DrawTriMeshNorm
(const std::vector<CVector3>& aP,
 const std::vector<int>& aTri)
{
  const int nTri = (int)aTri.size()/3;
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<nTri;itri++){
    const int i0 = aTri[itri*3+0];
    const int i1 = aTri[itri*3+1];
    const int i2 = aTri[itri*3+2];
    const CVector3& v0 = aP[i0];
    const CVector3& v1 = aP[i1];
    const CVector3& v2 = aP[i2];
    const CVector3& n = Normal(v0, v1, v2).Normalize();
    myGlNormal(n);
    myGlVertex(v0);
    myGlVertex(v1);
    myGlVertex(v2);
  }
  ::glEnd();
}

void DrawMeshTri_Edge
(const std::vector<CVector3>& aP,
 const std::vector<unsigned int>& aTri)
{
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ////
  /////
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  const int nTri = (int)aTri.size()/3;
  for (int itri = 0; itri<nTri; ++itri){
    const int i0 = aTri[itri*3+0];
    const int i1 = aTri[itri*3+1];
    const int i2 = aTri[itri*3+2];
    if( i0 == -1 ){
      assert(i1==-1); assert(i2==-1);
      continue;
    }
    assert(i0>=0&&i0<(int)aP.size());
    assert(i1>=0&&i1<(int)aP.size());
    assert(i2>=0&&i2<(int)aP.size());
    ::myGlVertex(aP[i0]); ::myGlVertex(aP[i1]);
    ::myGlVertex(aP[i1]); ::myGlVertex(aP[i2]);
    ::myGlVertex(aP[i2]); ::myGlVertex(aP[i0]);
  }
  ::glEnd();
  
  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}


void DrawMeshQuad_Face
(const std::vector<CVector3>& aPoint,
 const std::vector<unsigned int>& aQuad)
{
  ::glBegin(GL_QUADS);
  for(int iq=0;iq<(int)aQuad.size()/4;++iq){
    int iv0 = aQuad[iq*4+0];
    int iv1 = aQuad[iq*4+1];
    int iv2 = aQuad[iq*4+2];
    int iv3 = aQuad[iq*4+3];
    {
      CVector3 v01 = aPoint[iv1]-aPoint[iv0];
      CVector3 v12 = aPoint[iv2]-aPoint[iv1];
      CVector3 n = (v01^v12).Normalize();
      myGlNormal(n);
    }
    myGlVertex(aPoint[iv0]);
    myGlVertex(aPoint[iv1]);
    myGlVertex(aPoint[iv2]);
    myGlVertex(aPoint[iv3]);
  }
  ::glEnd();
}

void DrawPoint3D
(const std::vector<CVector3>& aPoint)
{
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_POINTS);
  for(int ip=0;ip<(int)aPoint.size();++ip){
    ::myGlVertex(aPoint[ip]);
  }
  ::glEnd();
}

void DrawQuad3D_Edge
(const std::vector<CVector3>& aPoint,
 const std::vector<unsigned int>& aQuad)
{
  ::glBegin(GL_LINES);
  for(int iq=0;iq<(int)aQuad.size()/4;++iq){
    int iv0 = aQuad[iq*4+0];
    int iv1 = aQuad[iq*4+1];
    int iv2 = aQuad[iq*4+2];
    int iv3 = aQuad[iq*4+3];
    myGlVertex(aPoint[iv0]);  myGlVertex(aPoint[iv1]);
    myGlVertex(aPoint[iv1]);  myGlVertex(aPoint[iv2]);
    myGlVertex(aPoint[iv2]);  myGlVertex(aPoint[iv3]);
    myGlVertex(aPoint[iv3]);  myGlVertex(aPoint[iv0]);
  }
  ::glEnd();
}

void DrawSingleHex_Edge
(const CVector3& p0, const CVector3& p1, const CVector3& p2, const CVector3& p3,
 const CVector3& p4, const CVector3& p5, const CVector3& p6, const CVector3& p7)
{
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
  myGlVertex(p0); myGlVertex(p4);
  myGlVertex(p1); myGlVertex(p5);
  myGlVertex(p2); myGlVertex(p6);
  myGlVertex(p3); myGlVertex(p7);
  ::glEnd();
}


void DrawGrid2D
(int ndivx, int ndivy,
 const CVector3& ex, const CVector3& ey, const CVector3& org)
{
  const CVector3& p00 = org;
  const CVector3& p10 = org + ex*ndivx;
  const CVector3& p01 = org + ey*ndivy;
  ::glBegin(GL_LINES);
  for(int ix=0;ix<ndivx+1;++ix){
    const CVector3& dx = ix*ex;
    myGlVertex(p00+dx);
    myGlVertex(p01+dx);
  }
  for(int iy=0;iy<ndivy+1;++iy){
    const CVector3& dy = iy*ey;
    myGlVertex(p00+dy);
    myGlVertex(p10+dy);
  }
  ::glEnd();
}

void DrawGridOutside
(int ndivx, int ndivy, int ndivz,
 double elen,
 const CVector3& org)
{
  DrawGrid2D(ndivx,ndivy, CVector3(elen,0,0), CVector3(0,elen,0), org);
  DrawGrid2D(ndivx,ndivy, CVector3(elen,0,0), CVector3(0,elen,0), org+CVector3(0,0,elen*ndivz));
  DrawGrid2D(ndivy,ndivz, CVector3(0,elen,0), CVector3(0,0,elen), org);
  DrawGrid2D(ndivy,ndivz, CVector3(0,elen,0), CVector3(0,0,elen), org+CVector3(elen*ndivx,0,0));
  DrawGrid2D(ndivz,ndivx, CVector3(0,0,elen), CVector3(elen,0,0), org);
  DrawGrid2D(ndivz,ndivx, CVector3(0,0,elen), CVector3(elen,0,0), org+CVector3(0,elen*ndivy,0));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isPickCircle
(const CVector2& sp,
 const CVector3& p,
 const CVector3& axis,
 double r,
 const float* mMV,
 const float* mPj,
 double pick_tol)
{
  const int ndiv = 32;
  double rdiv = 3.1415*2.0/ndiv;
  CVector3 x,y; GetVertical2Vector(axis, x, y);
  for(int idiv=0;idiv<ndiv+1;idiv++){
    int jdiv = idiv+1;
    CVector3 p0 = p+(r*sin(rdiv*idiv))*x+(r*cos(rdiv*idiv))*y;
    CVector3 p1 = p+(r*sin(rdiv*jdiv))*x+(r*cos(rdiv*jdiv))*y;
    CVector2 sp0 = screenXYProjection(p0, mMV, mPj);
    CVector2 sp1 = screenXYProjection(p1, mMV, mPj);
    double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
    if( sdist < pick_tol ){ return true; }
  }
  return false;
}

bool isPickCircle
(const CVector3& axis,
 const CVector3& org,
 double rad,
 const CVector3& src,
 const CVector3& dir,
 double pick_tol)
{
  double t = ((org-src)*axis)/(dir*axis);
  CVector3 p0 = src+t*dir;
  double rad0 = (p0-org).Length();
  if( fabs(rad-rad0) < pick_tol ) return true;
  return false;
}

double DragCircle
(const CVector2& sp0,
 const CVector2& sp1,
 const CVector3& p,
 const CVector3& axis,
 const float* mMV,
 const float* mPj)
{
  CVector2 spo0 = screenXYProjection(p, mMV, mPj);
  double area = TriArea(sp0, spo0, sp1);
  double angl = area / ( (sp0-spo0).Length() * (sp1-spo0).Length() );
  {
    CVector3 a3 = screenUnProjectionDirection(axis,mMV,mPj);
    if( a3.z < 0 ){ angl *= -1; }
  }
  return angl;
//  CMatrix3 R; R.SetRotMatrix_Cartesian(angl*axis);
//  return R;
}

bool isPickQuad
(const CVector3& p0,const CVector3& p1,const CVector3& p2,const CVector3& p3,
 const CVector2& sp, const CVector3& pick_dir,
 const float mMV[16], const float mPj[16],
 double eps)
{
  const CVector2 sp0 = screenXYProjection(p0, mMV, mPj);
  const CVector2 sp1 = screenXYProjection(p1, mMV, mPj);
  const CVector2 sp2 = screenXYProjection(p2, mMV, mPj);
  const CVector2 sp3 = screenXYProjection(p3, mMV, mPj);
  double a01 = TriArea(sp,sp0,sp1);
  double a12 = TriArea(sp,sp1,sp2);
  double a23 = TriArea(sp,sp2,sp3);
  double a30 = TriArea(sp,sp3,sp0);
  double a0123 = a01+a12+a23+a30;
  if( fabs(a0123) < 1.0e-10 ) return false;
  a01 /= a0123;
  a12 /= a0123;
  a23 /= a0123;
  a30 /= a0123;
  if( a01<eps || a12<eps || a23<eps || a30<eps ){ return false; }
  CVector3 n0123 = Normal(p0,p1,p2) + Normal(p1,p2,p3) + Normal(p2,p3,p0) + Normal(p3,p0,p1);
  if( n0123*pick_dir>0 ){ return false; }
  return true;
}

void DrawAxisHandler(double s, const CVector3& p)
{
  GLboolean is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1, 0, 0);
  ::DrawArrow(p,CVector3(+s, 0, 0));
  ::DrawArrow(p,CVector3(-s, 0, 0));
  
  ::glColor3d(0, 1, 0);
  ::DrawArrow(p, CVector3(0, +s, 0));
  ::DrawArrow(p, CVector3(0, -s, 0));
  
  ::glColor3d(0, 0, 1);
  ::DrawArrow(p, CVector3(0, 0, +s));
  ::DrawArrow(p, CVector3(0, 0, -s));
  
  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}

void DrawHandlerRotation
(const CVector3& pos,
 const double quat[4],
 double size,
 int ielem_picked)
{
  ::glDisable(GL_LIGHTING);
  if( ielem_picked == 0 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(1,0,0); }
  ::DrawCircleWire(QuatVec(quat,CVector3(1,0,0)), pos, size);
  if( ielem_picked == 1 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(0,1,0); }
  ::DrawCircleWire(QuatVec(quat,CVector3(0,1,0)), pos, size);
  if( ielem_picked == 2 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(0,0,1); }
  ::DrawCircleWire(QuatVec(quat,CVector3(0,0,1)), pos, size);
}

int PickHandlerRotation
(const CVector3& src, const CVector3& dir,
 const CVector3& pos, const double quat[4], double rad,
 double tol)
{
  CVector3 px,qx; Nearest_Line_Circle(px,qx, src,dir, pos,QuatVec(quat,CVector3(1,0,0)), rad);
  CVector3 py,qy; Nearest_Line_Circle(py,qy, src,dir, pos,QuatVec(quat,CVector3(0,1,0)), rad);
  CVector3 pz,qz; Nearest_Line_Circle(pz,qz, src,dir, pos,QuatVec(quat,CVector3(0,0,1)), rad);
  double dx = (px-src)*dir;
  double dy = (py-src)*dir;
  double dz = (pz-src)*dir;
  double lx = (px-qx).Length();
  double ly = (py-qy).Length();
  double lz = (pz-qz).Length();
  double dm = (fabs(dx)+fabs(dy)+fabs(dz))*1000;
  std::cout << lx << " " << ly << " " << lz << " " << dm << std::endl;
  if( lx>tol ){ dx = dm; }
  if( ly>tol ){ dy = dm; }
  if( lz>tol ){ dz = dm; }
  if( dx < dy && dx < dz  && dx < 0.9*dm ){ return 0; }
  if( dy < dz && dy < dx  && dy < 0.9*dm ){ return 1; }
  if( dz < dx && dz < dy  && dz < 0.9*dm ){ return 2; }
  return -1;
}


bool DragHandlerRot
(double quat[4], int ielem,
 const CVector2& sp0, const CVector2& sp1,
 const CVector3& pos,
 const float mMV[16], const float mPj[16])
{
  if( ielem>=0 && ielem<3 ){
    double vi[3] = {0,0,0}; vi[ielem] = 1;
    double vo[3]; QuatVec(vo, quat, vi);
    CVector3 v0(0,0,0); v0[ielem] = 1;
    CVector3 v1(vo[0],vo[1],vo[2]); v1.SetNormalizedVector();
    double ar = -DragCircle(sp0,sp1, pos, v1, mMV, mPj);
    double dq[4] = { cos(ar*0.5), v0.x*sin(ar*0.5), v0.y*sin(ar*0.5), v0.z*sin(ar*0.5) };
    double qtmp[4]; QuatMult(qtmp, dq, quat);
    QuatCopy(quat,qtmp);
    return true;
  }
  return false;
}

bool isPick_AxisHandler
(const CVector2& sp,
 const CVector3& p,
 const CVector3& axis,
 double len,
 const float* mMV,
 const float* mPj,
 double pick_tol)
{
  CVector2 sp0 = screenXYProjection(p+len*axis, mMV, mPj);
  CVector2 sp1 = screenXYProjection(p-len*axis, mMV, mPj);
  double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
  if (sdist < pick_tol){ return true; }
  return false;
}

CVector3 drag_AxisHandler
(const CVector2& sp0,
 const CVector2& sp1,
 const CVector3& p,
 const CVector3& axis,
 double len,
 const float* mMV,
 const float* mPj)
{
  CVector2 spa0 = screenXYProjection(p+len*axis, mMV, mPj);
  CVector2 spa1 = screenXYProjection(p-len*axis, mMV, mPj);
  double r = (spa0-spa1)*(sp1-sp0)/(spa0-spa1).SqLength();
  return r*axis*len;
}

bool isPickPoint
(const CVector2& sp,
 const CVector3& p,
 const float* mMV,
 const float* mPj,
 double pick_tol)
{
  CVector2 sp0 = screenXYProjection(p, mMV, mPj);
  if ((sp-sp0).Length() < pick_tol){ return true; }
  return false;
}


/////////////////////////////////

CVector3 QuatVec(const double quat[4], const CVector3& v0){
  const double v0a[3] = {v0.x,v0.y,v0.z};
  double v1a[3]; QuatVec(v1a,quat,v0a);
  return CVector3(v1a[0],v1a[1],v1a[2]);
}

CVector3 QuatConjVec(const double quat[4], const CVector3& v0){
  const double v0a[3] = {v0.x,v0.y,v0.z};
  double v1a[3]; QuatConjVec(v1a,quat,v0a);
  return CVector3(v1a[0],v1a[1],v1a[2]);
}


