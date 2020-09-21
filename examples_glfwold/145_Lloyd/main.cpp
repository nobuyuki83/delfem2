/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cstdlib>
#include <set>
#include <random>
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/color.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"

namespace dfm2 = delfem2;

// -----------------------------

void Step(
    std::vector<double>& aXY,
    const unsigned int ndiv,
    const std::vector<double>& aD)
{
  assert( aD.size() == ndiv*ndiv );
  const unsigned int np = aXY.size()/2;
  std::vector<unsigned int> aV; // Volonoi, index of point, distance for every pixels
  aV.resize(ndiv*ndiv );
  for(unsigned int ih=0; ih < ndiv; ++ih) {
    for(unsigned int iw=0;iw<ndiv;++iw){
      const double x0 = (iw + 0.5) / ndiv;
      const double y0 = 1.0 - (ih + 0.5) / ndiv;
      double min_dist = 4.0;
      unsigned int min_ip = UINT_MAX;
      for(unsigned int ip=0;ip<np;++ip){
        const double x1 = aXY[ip*2+0];
        const double y1 = aXY[ip*2+1];
        const double d01 = (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1);
        if( d01 > min_dist ){ continue; }
        min_dist = d01;
        min_ip = ip;
      }
      aV[ih*ndiv+iw] = min_ip;
    }
  }

  std::vector<double> awpw(np*3,1.0e-10);
  for(unsigned int ih=0; ih < ndiv; ++ih) {
    for (unsigned int iw = 0; iw < ndiv; ++iw) {
      const double x0 = (iw + 0.5) / ndiv;
      const double y0 = 1.0 - (ih + 0.5) / ndiv;
      double w0 = aD[ih*ndiv+iw];
      const unsigned int ip0 = aV[ih*ndiv+iw];
      awpw[ip0*3+0] += w0*x0;
      awpw[ip0*3+1] += w0*y0;
      awpw[ip0*3+2] += w0;
    }
  }

  for(unsigned int ip=0;ip<np;++ip){
    aXY[ip*2+0] = awpw[ip*3+0] / awpw[ip*3+2];
    aXY[ip*2+1] = awpw[ip*3+1] / awpw[ip*3+2];
  }
}

int main(int argc,char* argv[])
{
  const unsigned int ndiv = 128;
  const unsigned int np = 100;
  std::vector<double> aXY;
  std::vector<double> aDensity;
  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_real_distribution<double> dist(0,1);

  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.nav.camera.trans[0] = -0.5;
  viewer.nav.camera.trans[1] = -0.5;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 0.7;
  while ( true )
  {
    for(unsigned int istep=0;istep<200;++istep) {
      if( istep % 100  == 0 ) {
        aXY.resize(np*2);
        for (unsigned int ip = 0; ip < np; ++ip) {
          aXY[ip*2+0] = dist(rdeng);
          aXY[ip*2+1] = dist(rdeng);
        }
      }
      if( istep == 0 ){
        aDensity.assign(ndiv*ndiv,1.0);
      }
      if( istep == 100 ){
        aDensity.resize(ndiv*ndiv);
        for(unsigned int idiv=0;idiv<ndiv;++idiv){
          for(unsigned int jdiv=0;jdiv<ndiv;++jdiv) {
            double w = sin((double)idiv/ndiv*M_PI*2)*cos((double)jdiv/ndiv*M_PI*2);
            aDensity[idiv*ndiv+jdiv] = w*w;
          }
        }
      }
      Step(aXY,
           ndiv, aDensity);
      viewer.DrawBegin_oldGL();
      ::glColor3d(0, 0, 0);
      ::glPointSize(5);
      dfm2::opengl::DrawPoints2d_Points(aXY);
      double p0[2] = {0,0}, p1[2] = {1,1};
      dfm2::opengl::DrawBox2_Edge(p0,p1);
      viewer.DrawEnd_oldGL();
    }
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
