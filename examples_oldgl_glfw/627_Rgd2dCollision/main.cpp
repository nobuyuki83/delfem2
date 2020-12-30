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
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ------------------------------------------------------

class CRigidState2
{
public:
  bool is_fix;
  dfm2::CVec2d posg;
  double theta;
  dfm2::CVec2d velo;
  double omega;
  std::vector<dfm2::CVec2d> shape;
  dfm2::CVec2d posl;
};

void Draw(
    const dfm2::CVec2d& pos,
    double theta,
    const std::vector<dfm2::CVec2d>& shape)
{
  double mR[9]; dfm2::Mat3_AffineRotation(mR,theta);
  double mT[9]; dfm2::Mat3_AffineTranslation(mT,pos.p);
  double mTR[9]; dfm2::MatMat3(mTR,mT,mR);
  //
  ::glBegin(GL_LINES);
  ::glColor3d(0,0,0);
  for(unsigned int i0=0;i0<shape.size();++i0){
    unsigned int i1 = (i0+1)%shape.size();
    double p0[2]; dfm2::Vec2_Mat3Vec2_AffineProjection(p0,mTR,shape[i0].p);
    double p1[2]; dfm2::Vec2_Mat3Vec2_AffineProjection(p1,mTR,shape[i1].p);
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

  CRigidState2 rs;
  {
    rs.is_fix = false;
    rs.shape = {
        {-0.2,-0.1},
        {+0.2,-0.1},
        {+0.2,+0.1},
        {-0.2,+0.1}};
//    dfm2::SecondMomentOfArea_Polygon(rs.posl,area,pa1,I1,pa2,I2,rs.shape);
    rs.omega = 0.0;
    rs.theta = M_PI / 16;
    rs.velo = dfm2::CVec2d(0, 0);
    rs.posg = dfm2::CVec2d(0, 0.5);
  }

  while(true){
    //
    viewer.DrawBegin_oldGL();
    Draw(rs.posg,rs.theta,rs.shape);
    ::glBegin(GL_LINES);
    ::glVertex2d(-1,-0.5);
    ::glVertex2d(+1,-0.5);
    ::glEnd();
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ break; }
  }
  viewer.ExitIfClosed();
  return 0;
}


