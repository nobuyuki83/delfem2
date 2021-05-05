/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/rgd_v2m3.h"
#include "delfem2/geoplygn2_v2.h"
#include "delfem2/vec2.h"
#include "delfem2/mat3.h"
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ------------------------------------------------------

void Draw(
    const delfem2::CRigidState2& rs)
{
  const dfm2::CMat3d mT1RT0 = delfem2::rgd_v2m3::Mat3_Affine(rs.posl,rs.theta,rs.posg);
  ::glBegin(GL_LINES);
  ::glColor3d(0,0,0);
  for(unsigned int i0=0;i0<rs.shape.size();++i0){
    unsigned int i1 = (i0+1)%rs.shape.size();
    double p0[2]; dfm2::Vec2_Mat3Vec2_AffineProjection(p0,mT1RT0.mat,rs.shape[i0].p);
    double p1[2]; dfm2::Vec2_Mat3Vec2_AffineProjection(p1,mT1RT0.mat,rs.shape[i1].p);
    ::glVertex2dv(p0);
    ::glVertex2dv(p1);
  }
  ::glEnd();
  //
  ::glPointSize(5);
  ::glBegin(GL_POINTS);
  ::glColor3d(0,0,0);
  for(unsigned int i0=0;i0<rs.shape.size();++i0) {
    double p0[2];
    dfm2::Vec2_Mat3Vec2_AffineProjection(p0, mT1RT0.mat, rs.shape[i0].p);
    ::glVertex2dv(p0);
  }
  ::glEnd();
}

int main(int argc,char* argv[])
{
  dfm2::glfw::InitGLOld();
  dfm2::glfw::CViewer3 viewer;
  viewer.InitGL();
  viewer.camera.view_height = 1.5;

  std::vector<dfm2::CRigidState2> aRS(2);
  {
    dfm2::CRigidState2& rs = aRS[0];
    rs.is_fix = true;
    rs.shape = {
        {-1.0, 0.0},
        {+1.0, 0.0},
        {+1.0, 0.8},
        {+0.9, 0.8},
        {+0.9, 0.2},
        {-0.9, 0.2},
        {-0.9, 0.8},
        {-1.0, 0.8},
    };
  }
  {
    dfm2::CRigidState2& rs = aRS[1];
    rs.is_fix = false;
    unsigned int ndiv = 10;
    rs.shape.reserve((ndiv+1)*2);
    for(int i=0;i<ndiv+1;++i){
      rs.shape.emplace_back(i*0.4/ndiv,0.0 );
    }
    for(int i=(int)ndiv;i>=0;--i){
      rs.shape.emplace_back(i*0.4/ndiv,0.2 );
    }
  }

  for(dfm2::CRigidState2& rs : aRS){
    dfm2::CgArea_Polygon(rs.posl,rs.mass, rs.shape);
    rs.I = dfm2::RotationalMomentPolar_Polygon2(rs.shape,rs.posl);
    double rho = 1.0;
    rs.mass *= rho;
    rs.I *= rho;
    rs.omega = 0.0;
    rs.theta = 0;
    rs.velo = dfm2::CVec2d(0, 0);
    rs.posg = dfm2::CVec2d(0, 0.0);
  }

  aRS[1].posg.p[1] = 0.5;
  aRS[1].theta = 0.0;

  const dfm2::CVec2d gravity(0,-10);
  double dt = 0.01;

  double time = 0.0;
  while(true){
    {
      auto& aP = aRS[1].shape;
      unsigned int np = aP.size()/2;
      for(int ip=0;ip<np;++ip) {
        double x0 = aP[ip].x();
        aP[ip].p[1] = 0.05 * sin(x0 * 20 + time * 5);
      }
      for(int ip=0;ip<np;++ip){
        double x0 = aP[ip+np].x();
        aP[ip+np].p[1] = 0.05*sin(x0*20+time*5)+0.2;
      }
    }
    Steptime_Rgd2(aRS, dt, gravity);
    time += dt;
    //
    viewer.DrawBegin_oldGL();
    for(const dfm2::CRigidState2& rs : aRS) {
      Draw(rs);
    }
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ break; }
  }
  viewer.ExitIfClosed();
  return 0;
}


