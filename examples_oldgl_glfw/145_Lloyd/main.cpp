/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/sampler.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <cstdlib>
#include <set>
#include <random>


#define STB_IMAGE_IMPLEMENTATION
#include "delfem2/../../3rd_party/stb_image.h"

namespace dfm2 = delfem2;

void MakeDensity(
  std::vector<double>& aDensity,
  unsigned int imode,
  const unsigned int ndiv)
{
  if( imode == 0 ){
    aDensity.assign(ndiv*ndiv,500.0);
  }
  else if( imode == 1 ){
    aDensity.resize(ndiv*ndiv);
    for(unsigned int idiv=0;idiv<ndiv;++idiv){
      for(unsigned int jdiv=0;jdiv<ndiv;++jdiv) {
        double w0 = sin((double)idiv/ndiv*M_PI*2)*cos((double)jdiv/ndiv*M_PI*2);
        aDensity[idiv*ndiv+jdiv] = 2000*w0*w0;
      }
    }
  }
  else if( imode == 2 ){
    int width, height;
    std::string name_img_in_test_inputs = "tesla.png";
    int channels;
    unsigned char *img = stbi_load((std::string(PATH_INPUT_DIR)+"/"+name_img_in_test_inputs).c_str(),
                                   &width, &height, &channels, 0);
    assert( width == 128 && height == 128 && channels == 3 );
    for(unsigned int ipix=0;ipix<ndiv*ndiv;++ipix){
      double a = (double)img[ipix*3+0] + (double)img[ipix*3+1] + (double)img[ipix*3+2];
      double d0 = 255.0*3 - a;
      aDensity[ipix] = d0*d0*d0*0.00005;
    }
    delete[] img;
  }
}

void Draw(
    const delfem2::opengl::CViewer_GLFW& viewer,
    const std::vector<double>& aXY,
    const double min_xy[2],
    const double max_xy[2])
{
  viewer.DrawBegin_oldGL();
  ::glColor3d(0, 0, 0);
  ::glPointSize(5);
  dfm2::opengl::DrawPoints2d_Points(aXY);
  dfm2::opengl::DrawBox2_Edge(min_xy,max_xy);
  viewer.SwapBuffers();
  glfwPollEvents();
}

int main(int argc,char* argv[])
{
  std::vector<double> aXY;
  const unsigned int ndiv = 128;
  std::vector<double> aDensity;
  const double min_xy[2] = {0,0};
  const double max_xy[2] = {1,1};
  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_real_distribution<double> dist_x(min_xy[0], max_xy[0]);
  std::uniform_real_distribution<double> dist_y(min_xy[1], max_xy[1]);
  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.camera.trans[0] = -0.5;
  viewer.camera.trans[1] = -0.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  viewer.Init_oldGL();
  viewer.camera.view_height = 0.7;
  while ( true )
  {
    for(unsigned int imode=0;imode<3;++imode){
      MakeDensity(aDensity,imode,ndiv);
      aXY.resize(0);
      int icnt_fail = 0;
      for(unsigned int idart=0;idart<1000;++idart){
        for(int itr=0;itr<10;++itr) {
          double x0 = dist_x(rdeng);
          double y0 = dist_x(rdeng);
          auto ix0 = (unsigned int) floor((x0 - min_xy[0]) / (max_xy[0] - min_xy[0]) * ndiv);
          auto iy0 = (unsigned int) floor((y0 - min_xy[1]) / (max_xy[1] - min_xy[1]) * ndiv);
          iy0 = ndiv - iy0;
          double dist01 = 1.0 / aDensity[iy0 * ndiv + ix0];
          bool flg_fail = false;
          for (unsigned int ip = 0; ip < aXY.size() / 2; ++ip) {
            const double x1 = aXY[ip * 2 + 0];
            const double y1 = aXY[ip * 2 + 1];
            const double d01 = (x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1);
            if (d01 * d01 < dist01 * dist01) {
              flg_fail = true;
              icnt_fail++;
              break;
            }
          }
          if (!flg_fail) {
            icnt_fail = 0;
            aXY.push_back(x0);
            aXY.push_back(y0);
          }
        }
        Draw(viewer,aXY,min_xy,max_xy);
        if( icnt_fail > 100 ){ break; }
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
      for(int itr=0;itr<50;++itr){
        dfm2::Step_Lloyd2(aXY,
                          ndiv, aDensity, min_xy, max_xy);
        Draw(viewer,aXY,min_xy,max_xy);
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
    }
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
