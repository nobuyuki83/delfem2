/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/vec2.h"
#include "delfem2/geoplygn2_v2.h"
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ------------------------------------------------------

class CRigidState2
{
public:
  bool is_fix;
  std::vector<dfm2::CVec2d> shape;
  // below: derived from shape
  dfm2::CVec2d posl;
  double mass;
  double I;
  //
  // below: change over time
  dfm2::CVec2d posg;
  double theta;
  dfm2::CVec2d velo;
  double omega;
  // below: temp data for PBD
  dfm2::CVec2d posg_tmp;
  double theta_tmp;
};

void Draw(
    const CRigidState2& rs)
{
  dfm2::CMat3d mT0; dfm2::Mat3_AffineTranslation(mT0.mat,(-rs.posl).p);
  dfm2::CMat3d mR; dfm2::Mat3_AffineRotation(mR.mat,rs.theta);
  dfm2::CMat3d mT1; dfm2::Mat3_AffineTranslation(mT1.mat,rs.posg.p);
  dfm2::CMat3d mT1RT0 = mT1*mR*mT0;
  //
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
}

int main(int argc,char* argv[])
{
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 1.5;

  std::vector<CRigidState2> aRS(3);
  {
    CRigidState2& rs = aRS[0];
    rs.is_fix = false;
    rs.shape = { {0.0, 0.0}, {0.4, 0.0}, {0.4, 0.2}, {0.0, 0.2} };
  }
  {
    CRigidState2& rs = aRS[1];
    rs.is_fix = true;
    rs.shape = { {-1.0, 0.0}, {+1.0, 0.0}, {+1.0, 0.2}, {-1.0, 0.2} };
  }
  {
    CRigidState2& rs = aRS[2];
    rs.is_fix = false;
    rs.shape = { {0.0, 0.0}, {0.4, 0.0}, {0.4, 0.2}, {0.0, 0.2} };
  }

  for(CRigidState2& rs : aRS){
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

  aRS[0].posg.p[1] = 0.5;
  aRS[0].theta = M_PI*0.1;

  aRS[2].posg.p[0] = 0.1;
  aRS[2].posg.p[1] = 1.0;

  const dfm2::CVec2d gravity(0,-10);
  double dt = 0.001;

  while(true){
    for(CRigidState2& rs : aRS){
      if( rs.is_fix ) { continue; }
      rs.velo += gravity * dt;
      rs.posg_tmp = rs.posg + dt * rs.velo + dt * dt * gravity;
      rs.theta_tmp = rs.theta + dt * rs.omega;
    }
    for(unsigned int ir=0;ir<aRS.size();++ir) {
      CRigidState2& ri = aRS[ir];
      for(dfm2::CVec2d& pi : ri.shape) {
        const dfm2::CVec2d qi = (pi-ri.posl).Rotate(ri.theta_tmp) + ri.posg_tmp;
        for(unsigned int jr=0;jr<aRS.size();++jr) {
          CRigidState2& rj = aRS[jr];
          if( ir == jr ){ continue; }
          const unsigned int nej = rj.shape.size();
          for(unsigned int iej=0;iej<nej;++iej){
            const dfm2::CVec2d qj0 = (rj.shape[(iej+0)%nej]-rj.posl).Rotate(rj.theta_tmp) + rj.posg_tmp;
            const dfm2::CVec2d qj1 = (rj.shape[(iej+1)%nej]-rj.posl).Rotate(rj.theta_tmp) + rj.posg_tmp;
            const dfm2::CVec2d qje = dfm2::GetNearest_LineSeg_Point(qi,qj0,qj1);
            if( (qi - qje).Length() > 1.0e-2 ){ continue; }
            const dfm2::CVec2d nj = -(qj1-qj0).Rotate(M_PI*0.5).Normalize();
            const double penetration = (qje-qi)*nj;
            if( penetration < 0 ){ continue; }
            double deno = 0.0;
            if( !ri.is_fix ){
              deno += 1/ri.mass;
              double t0 = (qi-ri.posg).Rotate(M_PI*0.5) * nj;
              deno += t0*t0/ri.I;
            }
            if( !rj.is_fix ){
              deno += 1/rj.mass;
              double t0 = (qje-rj.posg).Rotate(M_PI*0.5) * (-nj);
              deno += t0*t0/rj.I;
            }
            deno = penetration/deno;
            if( !ri.is_fix ){
              ri.posg_tmp += (deno/ri.mass)*nj;
              ri.theta_tmp += deno * (qi-ri.posg).Rotate(M_PI*0.5) * nj / ri.I;
            }
            if( !rj.is_fix ){
              rj.posg_tmp += (deno/rj.mass)*(-nj);
              rj.theta_tmp += deno * (qje-rj.posg).Rotate(M_PI*0.5) * (-nj) / rj.I;
            }
          }
        }
      }
    }
    for(CRigidState2& rs : aRS) {
      if( rs.is_fix ) { continue; }
      rs.velo = (rs.posg_tmp - rs.posg)/dt;
      rs.omega = (rs.theta_tmp - rs.theta)/dt;
      rs.posg = rs.posg_tmp;
      rs.theta = rs.theta_tmp;
    }
    //
    viewer.DrawBegin_oldGL();
    for(const CRigidState2& rs : aRS) {
      Draw(rs);
    }
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ break; }
  }
  viewer.ExitIfClosed();
  return 0;
}


