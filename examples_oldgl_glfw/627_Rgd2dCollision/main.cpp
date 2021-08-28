/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/rgd_v2m3.h"
#include "delfem2/geoplygn2_v2.h"
#include "delfem2/vec2.h"
#include "delfem2/mat3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// -----------------------------------------------------

void Draw(
    const delfem2::CRigidState2 &rs) {
  const dfm2::CMat3d mT1RT0 = delfem2::rgd_v2m3::Mat3_Affine(rs.posl, rs.theta, rs.posg);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (unsigned int i0 = 0; i0 < rs.shape.size(); ++i0) {
    unsigned int i1 = (i0 + 1) % rs.shape.size();
    double p0[2];
    dfm2::Vec2_Mat3Vec2_AffineProjection(p0, mT1RT0.mat, rs.shape[i0].p);
    double p1[2];
    dfm2::Vec2_Mat3Vec2_AffineProjection(p1, mT1RT0.mat, rs.shape[i1].p);
    ::glVertex2dv(p0);
    ::glVertex2dv(p1);
  }
  ::glEnd();
}

int main(
    [[maybe_unused]] int argc,
    [[maybe_unused]] char *argv[]) {
  dfm2::glfw::InitGLOld();
  dfm2::glfw::CViewer3 viewer;
  viewer.InitGL();
  viewer.camera.view_height = 1.5;

  std::vector<dfm2::CRigidState2> aRS(4);
  {
    dfm2::CRigidState2 &rs = aRS[0];
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
    dfm2::CRigidState2 &rs = aRS[1];
    rs.is_fix = false;
    rs.shape = {{0.0, 0.0}, {0.4, 0.0}, {0.4, 0.2}, {0.0, 0.2}};
  }
  {
    dfm2::CRigidState2 &rs = aRS[2];
    rs.is_fix = false;
    rs.shape = {{0.0, 0.0}, {0.4, 0.0}, {0.4, 0.2}, {0.0, 0.2}};
  }
  {
    dfm2::CRigidState2 &rs = aRS[3];
    rs.is_fix = false;
    rs.shape = {{0.0, 0.0}, {0.4, 0.0}, {0.4, 0.2}, {0.0, 0.2}};
  }

  for (dfm2::CRigidState2 &rs : aRS) {
    dfm2::CgArea_Polygon(rs.posl, rs.mass, rs.shape);
    rs.I = dfm2::RotationalMomentPolar_Polygon2(rs.shape, rs.posl);
    double rho = 1.0;
    rs.mass *= rho;
    rs.I *= rho;
    rs.omega = 0.0;
    rs.theta = 0;
    rs.velo = dfm2::CVec2d(0, 0);
    rs.posg = dfm2::CVec2d(0, 0.0);
    rs.shape_velo.resize(rs.shape.size(), dfm2::CVec2d(0, 0));
  }

  aRS[1].posg.p[1] = 0.5;
  aRS[1].theta = M_PI * 0.1;

  aRS[2].posg.p[0] = 0.1;
  aRS[2].posg.p[1] = 1.0;

  aRS[3].posg.p[0] = -0.2;
  aRS[3].posg.p[1] = 1.3;

  const dfm2::CVec2d gravity(0, -10);
  double dt = 0.005;

  while (true) {
    Steptime_Rgd2(aRS, dt, gravity);
    //
    viewer.DrawBegin_oldGL();
    for (const dfm2::CRigidState2 &rs : aRS) {
      Draw(rs);
    }
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  viewer.ExitIfClosed();
  return 0;
}


