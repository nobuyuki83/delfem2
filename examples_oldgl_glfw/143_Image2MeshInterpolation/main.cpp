/*
 * Copyright (c) 2020-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <vector>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/img_ioppm.h"
#include "delfem2/img_interpolation.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

namespace dfm2 = delfem2;

void MeshTri2D_Square(
    std::vector<double>& aXY,
    std::vector<unsigned int>& aTri,
    double elen)
{
  namespace dfm2 = delfem2;
  //
  std::vector< std::vector<double> > aaXY0;
  {
    std::vector<double> aXY0;
    aXY0.push_back(0); aXY0.push_back(0);
    aXY0.push_back(1); aXY0.push_back(0);
    aXY0.push_back(1); aXY0.push_back(1);
    aXY0.push_back(0); aXY0.push_back(1);
    aaXY0.push_back(aXY0);
  }

  std::vector<dfm2::CDynPntSur> aPo2D;
  std::vector<dfm2::CVec2d> aVec2;
  std::vector<dfm2::CDynTri> aETri;
  { // initial tesselation
    std::vector<int> loopIP_ind, loopIP;
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY0);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if( elen > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     elen );
    }
    Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                   loopIP_ind,loopIP);
  }
  AssertDTri(aETri);
  AssertMeshDTri(aPo2D,aETri);
  CheckTri(aPo2D,aETri,aVec2);

  // filling inside
  if( elen > 1.0e-10 ){
    dfm2::CInputTriangulation_Uniform param(1.0);
    std::vector<int> aFlgPnt(aPo2D.size());
    std::vector<unsigned int> aFlgTri(aETri.size(),0);
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                  aVec2.size(), 0, elen, param);
  }
  AssertDTri(aETri);
  AssertMeshDTri(aPo2D,aETri);
  CheckTri(aPo2D,aETri,aVec2);

  // export to stl array
  dfm2::MeshTri2D_Export(aXY,aTri, aVec2, aETri);
}


// --------------------------------------

int main()
{
  std::vector<double> vec_xy;
  std::vector<unsigned int> vec_idx_tri;
  MeshTri2D_Square(
      vec_xy, vec_idx_tri,
      0.02);
  std::cout << "npoint: " << vec_xy.size()/2 << std::endl;

  std::vector<double> aColor;
  {
    aColor.resize(vec_xy.size()/2*3, 0.0);
    int width, height;
    std::string name_img_in_test_inputs = "lenna.png";
    // -----
    int channels;
    unsigned char *img = stbi_load(
        (std::string(PATH_INPUT_DIR)+"/"+name_img_in_test_inputs).c_str(),
        &width, &height, &channels, 0);
    std::cout << width << " " << height << " " << channels << std::endl;
    unsigned int nXY = static_cast<unsigned int>(vec_xy.size()/2);
    aColor.resize(nXY*3);
    for(unsigned int ip=0;ip<nXY;++ip) {
      double x = vec_xy[ip * 2 + 0] * (width - 1);
      double y = (1.0 - vec_xy[ip * 2 + 1]) * (height - 1);
      dfm2::ImageInterpolation_Bilinear(
          aColor.data() + ip * 3,
          width, height, img,
          x, y);
    }
    stbi_image_free(img);
  }
  std::cout << vec_xy.size()/2 << " " << 220*220 << std::endl;


  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.trans[0] = -0.5;
  viewer.camera.trans[1] = -0.5;
  while (!glfwWindowShouldClose(viewer.window))
  {
    // ---------
    viewer.DrawBegin_oldGL();
    ::glColor3d(0,0,0);
    dfm2::opengl::DrawMeshTri2D_Edge(vec_idx_tri, vec_xy);
//    ::glColor3d(1,1,1);
//    dfm2::opengl::DrawMeshTri2D_Face(aTri,aXY);
    dfm2::opengl::DrawMeshTri2D_FaceColor(
        vec_idx_tri.data(), static_cast<unsigned int>(vec_idx_tri.size()/3),
        vec_xy.data(), static_cast<unsigned int>(vec_xy.size()/2),
        aColor.data());
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
