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
#include "delfem2/funcs_gl.h"
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
  this->nResX = nresw;
  this->nResY = nresh;
  this->lengrid = elen;
  this->z_range = depth_max;
  z_axis[0] = dir_prj[0];  z_axis[1] = dir_prj[1];  z_axis[2] = dir_prj[2];
  origin[0] = org_prj[0];  origin[1] = org_prj[1];  origin[2] = org_prj[2];
  x_axis[0] = dir_width[0];  x_axis[1] = dir_width[1];  x_axis[2] = dir_width[2];
}

void CDepth::SetView(){
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  ViewTransformation(x_axis,z_axis,origin);
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::glOrtho(0.0, +lengrid*nResX,
            0.0, +lengrid*nResY,
            0, z_range);
  ::glMatrixMode(GL_MODELVIEW);
}

void CDepth::Start()
{
  glGetIntegerv(GL_VIEWPORT, view); // current viewport
  ::glViewport(0, 0, nResX, nResY);
  
  ::glClearColor(0.0f, 0.0f, 0.0f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glDisable(GL_BLEND);
  ::glEnable(GL_DEPTH_TEST);
  this->SetView();
}

void CDepth::End()
{
  //  glFinish();
  aZ.resize(nResX*nResY);
  
  glReadBuffer(GL_DEPTH_ATTACHMENT);
  glReadPixels(0, 0, nResX, nResY, GL_DEPTH_COMPONENT, GL_FLOAT, aZ.data());
  int n = nResX*nResY;
  for(int i=0;i<n;++i){ aZ[i] *= (-1.0*z_range); }
  
  ::glViewport(view[0], view[1], view[2], view[3]);
}

void CDepth::Draw() const {
  
  ::glPointSize(1);
  this->Draw_Point();
  /////
  ::glLineWidth(3);
  this->Draw_Axis();
  ////
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  this->Draw_BoundingBox();
}

void CDepth::Draw_Axis() const
{
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ModelTransformation(x_axis, z_axis, origin);
  DrawAxis(draw_len_axis);
    ::glPopMatrix();
}

void CDepth::Draw_BoundingBox() const
{
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ModelTransformation(x_axis, z_axis, origin);
  ::glLineWidth(3);
  Draw_AABB3D_MinMaxXYZ_Edge(0.0, lengrid*nResX, 0.0, lengrid*nResY, 0.0, -z_range);
  ::glPopMatrix();
}

void CDepth::Draw_Point() const
{
  ::glDisable(GL_LIGHTING);
  if( (int)aZ.size() != nResX*nResY ) return;
  if( color.size() == 3 ){ ::glColor3dv(color.data()); }
  if( color.size() == 4 ){ ::glColor4dv(color.data()); }
  ::glBegin(GL_POINTS);
  double o[3];
  for(int iy=0;iy<nResY;++iy){
    for(int ix=0;ix<nResX;++ix){
      this->getGPos(o,ix,iy,aZ[iy*nResX+ix]);
      ::glVertex3dv(o);
    }
  }
  ::glEnd();
}

void CDepth::getGPos(double p[3], int ix, int iy, double depth) const
{
  const CVector3& dx = x_axis;
  const CVector3& dz = z_axis;
  const CVector3& dy = Cross(dz,dx);
  double lz = depth;
  double lx = (ix+0.5)*lengrid;
  double ly = (iy+0.5)*lengrid;
  CVector3 vp = lx*dx+ly*dy+lz*dz + origin;
  p[0] = vp.x;
  p[1] = vp.y;
  p[2] = vp.z;
}
