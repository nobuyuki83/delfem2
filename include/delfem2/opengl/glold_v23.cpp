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


#include "delfem2/vec2.h"
#include "delfem2/vec3.h"

#include "delfem2/opengl/glold_v23.h"

namespace dfm2 = delfem2;

//----------------------------------------------------

void delfem2::opengl::myGlVertex(const CVector3& v)
{
  ::glVertex3d(v.x,v.y,v.z);
}

void delfem2::opengl::myGlTranslate(const CVector3& v)
{
  ::glTranslated(v.x,v.y,v.z);
}

void delfem2::opengl::myGlNormal(const CVector3& n)
{
  ::glNormal3d(n.x,n.y,n.z);
}

void delfem2::opengl::myGlNormal(const CVector3& a, const CVector3& b, const CVector3& c)
{
  CVector3 n; UnitNormal(n, a, b, c);
  ::glNormal3d(n.x,n.y,n.z);
}

void delfem2::opengl::myGlVertex(int i, const std::vector<CVector3>& aV)
{
  const CVector3& v = aV[i];
  opengl::myGlVertex(v);
}

void delfem2::opengl::myGlVertex2(int i, const std::vector<double>& vec)
{
  ::glVertex3d(vec[i*2], vec[i*2+1], +0.0);
}

void delfem2::opengl::myGlVertex3(
    unsigned int i,
    const std::vector<double>& vec)
{
  ::glVertex3d(vec[i*3+0], vec[i*3+1], vec[i*3+2]);
}

void delfem2::opengl::ModelTransformation
 (const CVector3& dx, const CVector3& dz, const CVector3& origin)
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

void delfem2::opengl::ViewTransformation
 (const CVector3& dx, const CVector3& dz, const CVector3& origin)
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

void delfem2::opengl::myGlVertex(
    unsigned int i,
    const std::vector<CVector2>& aP)
{
  ::glVertex3d(aP[i].x, aP[i].y, +0.0);
}

void delfem2::opengl::myGlVertex(const CVector2& v)
{
  ::glVertex2d(v.x,v.y);
}

//--------------------------------------------------------

void delfem2::opengl::DrawCylinderWire
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

void delfem2::opengl::DrawCylinder
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
      myGlNormal(n);
      myGlVertex(p0+r*sin((idiv+0)*dt)*x+r*cos((idiv+0)*dt)*y);
      myGlVertex(p1+r*sin((idiv+0)*dt)*x+r*cos((idiv+0)*dt)*y);
      myGlVertex(p1+r*sin((idiv+1)*dt)*x+r*cos((idiv+1)*dt)*y);
      myGlVertex(p0+r*sin((idiv+1)*dt)*x+r*cos((idiv+1)*dt)*y);
    }
    ::glEnd();
  }
  {
    ::glBegin(GL_TRIANGLES);
    for (int idiv = 0; idiv<ndivt; idiv++){
      CVector3 v0 = p1+(r*sin((idiv+0)*dt))*x+(r*cos((idiv+0)*dt))*y;
      CVector3 v1 = p1+(r*sin((idiv+1)*dt))*x+(r*cos((idiv+1)*dt))*y;
      const CVector3& v2 = p1;
      CVector3 n; UnitNormal(n, v1, v0, v2);
      myGlNormal(n);
      myGlVertex(v0);
      myGlVertex(v2);
      myGlVertex(v1);
    }
    ::glEnd();
  }
}


void delfem2::opengl::DrawArrow
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
      myGlNormal(n);
      myGlVertex(p0      +r0*sin((idiv+0)*dt)*x+r0*cos((idiv+0)*dt)*y);
      myGlVertex(p0+d*0.8+r0*sin((idiv+0)*dt)*x+r0*cos((idiv+0)*dt)*y);
      myGlVertex(p0+d*0.8+r0*sin((idiv+1)*dt)*x+r0*cos((idiv+1)*dt)*y);
      myGlVertex(p0      +r0*sin((idiv+1)*dt)*x+r0*cos((idiv+1)*dt)*y);
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
      myGlNormal(n);
      myGlVertex(v0);
      myGlVertex(v2);
      myGlVertex(v1);
    }
    ::glEnd();
  }
}

void delfem2::opengl::DrawCircleArrow
(const CVector3& org, CVector3 axis, double offset)
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
        myGlVertex(q0 + r0*sin((jdiv+0)*dt)*s0 + r0*cos((jdiv+0)*dt)*z);
        myGlVertex(q0 + r0*sin((jdiv+1)*dt)*s0 + r0*cos((jdiv+1)*dt)*z);
        myGlVertex(q1 + r0*sin((jdiv+1)*dt)*s1 + r0*cos((jdiv+1)*dt)*z);
        myGlVertex(q1 + r0*sin((jdiv+0)*dt)*s1 + r0*cos((jdiv+0)*dt)*z);
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
      const CVector3& v2 = q1;
      CVector3 n; UnitNormal(n, v0, v2, v1);
      ::glNormal3d(n.x,n.y,n.z);
      myGlVertex(v0);
      myGlVertex(v2);
      myGlVertex(v1);
    }
    ::glEnd();
  }
}

//--------------------------------------------------------

void delfem2::opengl::DrawCircleWire
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

void delfem2::opengl::DrawCircleSolid
(const CVector3& axis,
 const CVector3& org,
 double r)
{
  const CVector3& z = axis;
  CVector3 x, y; GetVertical2Vector(z, x, y);
  const int ndivt = 32;
  const double dt = 3.1415*2.0/ndivt;
  {
    ::glBegin(GL_TRIANGLES);
    for (int idiv = 0; idiv<ndivt; idiv++){
      CVector3 v0 = org+(r*sin((idiv+0)*dt))*x+(r*cos((idiv+0)*dt))*y;
      CVector3 v1 = org+(r*sin((idiv+1)*dt))*x+(r*cos((idiv+1)*dt))*y;
      const CVector3& v2 = org;
      CVector3 n; UnitNormal(n, v1, v0, v2);
      myGlNormal(n);
      myGlVertex(v0);
      myGlVertex(v2);
      myGlVertex(v1);
    }
    ::glEnd();
  }
}

void delfem2::opengl::DrawArcSolid
(const CVector3& axis,
 const CVector3& org,
 double ru, // rin
 double rv, // rout
 double rads,
 double rade)
{
  const CVector3& z = axis;
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
      myGlNormal(n);
      myGlVertex(v0);
      myGlVertex(v1);
      myGlVertex(u1);
      myGlVertex(u0);
    }
    ::glEnd();
  }
}

void delfem2::opengl::DrawSingleQuad_Edge
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

void delfem2::opengl::DrawSingleQuad_FaceNorm
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


void delfem2::opengl::drawPolyLine
(const std::vector<CVector2>& aP)
{
  ::glBegin(GL_LINES);
  for (size_t ip = 0; ip<aP.size()-1; ip++){
    unsigned int jp = ip+1;
    opengl::myGlVertex(ip,aP);
    opengl::myGlVertex(jp,aP);
  }
  ::glEnd();
  //
  ::glBegin(GL_POINTS);
  for (size_t ip = 0; ip<aP.size(); ip++){
    opengl::myGlVertex(ip,aP);
  }
  ::glEnd();
}

void delfem2::opengl::drawPolyLine3D
 (const std::vector<CVector3>& aP)
{
  if( aP.empty() ) return;
  ::glBegin(GL_LINES);
  for (size_t ip = 0; ip<aP.size()-1; ip++){
    unsigned int jp = ip+1;
    myGlVertex(ip,aP);
    myGlVertex(jp,aP);
  }
  ::glEnd();
  //
  ::glBegin(GL_POINTS);
  for (size_t ip = 0; ip<aP.size(); ip++){
    myGlVertex(ip,aP);
  }
  ::glEnd();
}

void delfem2::opengl::drawPolyLine2D
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

void delfem2::opengl::Draw_MeshTri
(const std::vector<CVector2>& aP,
 const std::vector<unsigned int>& aTri)
{
  const int nTri = (int)aTri.size()/3;
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<nTri;itri++){
    const int i0 = aTri[itri*3+0];
    const int i1 = aTri[itri*3+1];
    const int i2 = aTri[itri*3+2];
    const CVector2& v0 = aP[i0];
    const CVector2& v1 = aP[i1];
    const CVector2& v2 = aP[i2];
    myGlVertex(v0);
    myGlVertex(v1);
    myGlVertex(v2);
  }
  ::glEnd();
}

void delfem2::opengl::Draw_MeshTri_Edge
(const std::vector<CVector2>& aP,
 const std::vector<unsigned int>& aTri)
{
  //  const unsigned int nxys = (int)aXY.size()/2;
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  const int nTri = (int)aTri.size()/3;
  for(int itri=0;itri<nTri;itri++){
    const unsigned int i0 = aTri[itri*3+0];
    const unsigned int i1 = aTri[itri*3+1];
    const unsigned int i2 = aTri[itri*3+2];
    const CVector2& v0 = aP[i0];
    const CVector2& v1 = aP[i1];
    const CVector2& v2 = aP[i2];
    myGlVertex(v0);  myGlVertex(v1);
    myGlVertex(v1);  myGlVertex(v2);
    myGlVertex(v2);  myGlVertex(v0);
  }
  ::glEnd();
}

void delfem2::opengl::DrawTriMeshNorm
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

void delfem2::opengl::DrawMeshTri_Edge
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
    myGlVertex(aP[i0]); myGlVertex(aP[i1]);
    myGlVertex(aP[i1]); myGlVertex(aP[i2]);
    myGlVertex(aP[i2]); myGlVertex(aP[i0]);
  }
  ::glEnd();
  
  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}


void delfem2::opengl::DrawMeshQuad_Face
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

void delfem2::opengl::DrawPoint3D
(const std::vector<CVector3>& aPoint)
{
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_POINTS);
  for(const auto & p : aPoint){
    myGlVertex(p);
  }
  ::glEnd();
}

void delfem2::opengl::DrawQuad3D_Edge
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

void delfem2::opengl::DrawSingleHex_Edge
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


void delfem2::opengl::DrawGrid2D
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

void delfem2::opengl::DrawGridOutside
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

// -----------------------------------------------------------------

void delfem2::opengl::DrawAxisHandler(double s, const CVector3& p)
{
  GLboolean is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1, 0, 0);
  opengl::DrawArrow(p,CVector3(+s, 0, 0));
  opengl::DrawArrow(p,CVector3(-s, 0, 0));
  
  ::glColor3d(0, 1, 0);
  opengl::DrawArrow(p, CVector3(0, +s, 0));
  opengl::DrawArrow(p, CVector3(0, -s, 0));
  
  ::glColor3d(0, 0, 1);
  opengl::DrawArrow(p, CVector3(0, 0, +s));
  opengl::DrawArrow(p, CVector3(0, 0, -s));
  
  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}

void delfem2::opengl::DrawHandlerRotation_PosQuat
(const CVector3& pos,
 const double quat[4],
 double size,
 int ielem_picked)
{
  ::glDisable(GL_LIGHTING);
  {
    if( ielem_picked == 0 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(1,0,0); }
    const CVector3& ax = QuatVec(quat,CVector3(1,0,0));
    opengl::DrawCircleWire(ax, pos, size);
  }
  {
    if( ielem_picked == 1 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(0,1,0); }
    const CVector3& ay = QuatVec(quat,CVector3(0,1,0));
    opengl::DrawCircleWire(ay, pos, size);
  }
  {
    if( ielem_picked == 2 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(0,0,1); }
    const CVector3& az = QuatVec(quat,CVector3(0,0,1));
    opengl::DrawCircleWire(az, pos, size);
  }
}

void delfem2::opengl::DrawHandlerRotation_Mat4
(const double Mat[16],
 double size,
 int ielem_picked)
{
  ::glDisable(GL_LIGHTING);
  {
    if( ielem_picked == 0 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(1,0,0); }
    const CVector3& ax = Mat4Vec(Mat,CVector3(1,0,0)).Normalize();
    const CVector3 pos(Mat[3],Mat[7],Mat[11]);
    opengl::DrawCircleWire(ax, pos, size);
  }
  {
    if( ielem_picked == 1 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(0,1,0); }
    const CVector3& ay = Mat4Vec(Mat,CVector3(0,1,0)).Normalize();
    const CVector3 pos(Mat[3],Mat[7],Mat[11]);
    opengl::DrawCircleWire(ay, pos, size);
  }
  {
    if( ielem_picked == 2 ){ ::glColor3d(1,1,0); }   else{ ::glColor3d(0,0,1); }
    const CVector3& az = Mat4Vec(Mat,CVector3(0,0,1)).Normalize();
    const CVector3 pos(Mat[3],Mat[7],Mat[11]);
    opengl::DrawCircleWire(az, pos, size);
  }
}


