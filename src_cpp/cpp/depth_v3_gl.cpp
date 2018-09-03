#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stack>
#include <cstring>
#include <cstdlib>

#if defined(__APPLE__) && defined(__MACH__) // Mac
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
  #include <GL/glu.h>
#elif defined(WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else
  #include <GL/glew.h>
  #include <GL/gl.h>
#endif
////
#include "delfem2/mat3.h"
#include "delfem2/v23m3q.h"

#include "delfem2/v23_gl.h"  // vec3, mat3

#include "delfem2/depth_v3_gl.h"


/////////////////////////////////

void CDepth::SetColor(double r, double g, double b){
  color[0] = r;
  color[1] = g;
  color[2] = b;
}

void CDepth::SetCoord
(int nresw, int nresh, double elen, double depth_max,
 const std::vector<double>& org_prj,
 const std::vector<double>& dir_prj,
 const std::vector<double>& dir_width)
{
//  std::cout << nresw << " " << nresh << " " << elen << std::endl;
//  std::cout << "dir_prj: " << dirPrj[0] << " " << dirPrj[1] << " " << dirPrj[2] << std::endl;
//  std::cout << "org_prj: " << orgPrj[0] << " " << orgPrj[1] << " " << orgPrj[2] << std::endl;
//  std::cout << "width_prj: " << dirWidth[0] << " " << dirWidth[1] << " " << dirWidth[2] << std::endl;
  this->nResW = nresw;
  this->nResH = nresh;
  this->lengrid = elen;
  this->depth_max = depth_max;
  dirPrj[0] = dir_prj[0];  dirPrj[1] = dir_prj[1];  dirPrj[2] = dir_prj[2];
  orgPrj[0] = org_prj[0];  orgPrj[1] = org_prj[1];  orgPrj[2] = org_prj[2];
  dirWidth[0] = dir_width[0];  dirWidth[1] = dir_width[1];  dirWidth[2] = dir_width[2];
}


void CDepth::SetView(){
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  {
    const CVector3& dx = dirWidth;
    const CVector3& dz = dirPrj;
    const CVector3& dy = Cross(dz,dx);
    CMatrix3 R = Mat3(dx,dy,dz);
    CVector3 o = R.Trans()*CVector3(orgPrj);
    double A[16];
    A[ 0] = dx.x;  A[ 1] = dy.x;  A[ 2] = dz.x;  A[ 3] = 0;
    A[ 4] = dx.y;  A[ 5] = dy.y;  A[ 6] = dz.y;  A[ 7] = 0;
    A[ 8] = dx.z;  A[ 9] = dy.z;  A[10] = dz.z;  A[11] = 0;
    A[12] = -o.x;  A[13] = -o.y;  A[14] = -o.z;  A[15] = 1;
    ::glMultMatrixd(A);
  }
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::glOrtho(-lengrid*nResW*0.5, +lengrid*nResW*0.5,
            -lengrid*nResH*0.5, +lengrid*nResH*0.5,
            0, depth_max);
  ::glMatrixMode(GL_MODELVIEW);
}

void CDepth::Start()
{
  glGetIntegerv(GL_VIEWPORT, view); // current viewport
  ::glViewport(0, 0, nResW, nResH);
  
  ::glClearColor(0.0f, 0.0f, 0.0f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glDisable(GL_BLEND);
  ::glEnable(GL_DEPTH_TEST);
  this->SetView();
}

void CDepth::End()
{
  //  glFinish();
  aDepth.resize(nResW*nResH);
  
  glReadBuffer(GL_DEPTH_ATTACHMENT);
  glReadPixels(0, 0, nResW, nResH, GL_DEPTH_COMPONENT, GL_FLOAT, aDepth.data());
  
  ::glViewport(view[0], view[1], view[2], view[3]);
}

void CDepth::Draw() const {
  this->Draw_Point();
  const CVector3& dx = dirWidth;
  const CVector3& dz = dirPrj;
  const CVector3 dy = Cross(dz,dx);
  double p0[3] = {
    orgPrj[0]+dirWidth[0]*lengrid*nResW*0.5,
    orgPrj[1]+dirWidth[1]*lengrid*nResW*0.5,
    orgPrj[2]+dirWidth[2]*lengrid*nResW*0.5 };
  double p1[3] = {
    orgPrj[0]+dy.x*lengrid*nResH*0.5,
    orgPrj[1]+dy.y*lengrid*nResH*0.5,
    orgPrj[2]+dy.z*lengrid*nResH*0.5 };
  double p2[3] = {
    orgPrj[0]+dirPrj[0]*lengrid*depth_max*0.1,
    orgPrj[1]+dirPrj[1]*lengrid*depth_max*0.1,
    orgPrj[2]+dirPrj[2]*lengrid*depth_max*0.1 };
  ////////////////////////////
  ::glDisable(GL_LIGHTING);
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  ::glVertex3dv(orgPrj);
  ::glEnd();
  ////
  ::glLineWidth(3);
  ::glColor3d(1,0,0);
  ::glBegin(GL_LINES);
  ::glVertex3dv(orgPrj);
  ::glVertex3dv(p0);
  ::glColor3d(0,1,0);
  ::glVertex3dv(orgPrj);
  ::glVertex3dv(p1);
  ::glColor3d(0,0,1);
  ::glVertex3dv(orgPrj);
  ::glVertex3dv(p2);
  ::glEnd();
}

void CDepth::Draw_Point() const
{
  ::glDisable(GL_LIGHTING);
  if( (int)aDepth.size() != nResW*nResH ) return;
  if( color.size() == 3 ){ ::glColor3dv(color.data()); }
  if( color.size() == 4 ){ ::glColor4dv(color.data()); }
  ::glPointSize(1);
  ::glBegin(GL_POINTS);
  double o[3];
  for(int i=0;i<nResW*nResH;++i){
    this->getGPos(o,i,aDepth[i]);
    ::glVertex3dv(o);
  }
  ::glEnd();
}

void CDepth::getGPos(double p[3], int i, double depth) const{
  const CVector3& dx = dirWidth;
  const CVector3& dz = dirPrj;
  const CVector3 dy = Cross(dz,dx);
  CMatrix3 R = Mat3(dx,dy,dz);
  int ih = i/nResW;
  int iw = i-ih*nResW;
  double lz = -depth*depth_max;
  double lx = (iw+0.5)*lengrid-0.5*nResW*lengrid;
  double ly = (ih+0.5)*lengrid-0.5*nResH*lengrid;
  CVector3 vp = R*CVector3(lx,ly,lz) + orgPrj;
  p[0] = vp.x;
  p[1] = vp.y;
  p[2] = vp.z;
}
