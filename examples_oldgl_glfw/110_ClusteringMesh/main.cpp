/*
 * Copyright (c) 2019-2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cstdlib>
#include <random>
#include <set>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/dijkstra.h"
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshuni.h"
#include "delfem2/color.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"

namespace dfm2 = delfem2;

// ---------------------------

int main(
	[[maybe_unused]] int argc,
	[[maybe_unused]] char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;

  delfem2::Read_Ply(
//      std::string(PATH_INPUT_DIR)+"/bunny_34k.ply",
      std::string(PATH_INPUT_DIR)+"/arm_16k.ply",
      aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<unsigned int> aTriSuTri;
  ElSuEl_MeshElem(aTriSuTri,
      aTri.data(), aTri.size()/3,
      delfem2::MESHELEM_TRI,
      aXYZ.size()/3);

  // ------

  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_int_distribution<unsigned int> ncluster_gen(1,100);

  // above: data preparation
  // -----------------------
  // below: view

  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  viewer.camera.view_height = 0.5;
  viewer.camera.camera_rot_mode  = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    const unsigned int ncluster = ncluster_gen(rdeng);
    std::vector< std::pair<int,dfm2::CColor> > aColor;
    for(unsigned int ic=0;ic<ncluster;++ic){
      dfm2::CColor c;
      c.setRandomVividColor();
      aColor.emplace_back(2,c);
    }
    std::vector<unsigned int> aFlgElm;
    dfm2::MeshClustering(
		aFlgElm,ncluster,aTriSuTri,
		static_cast<unsigned int>(aTri.size()/3));
    //
    for(unsigned int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0, 0, 0);
      delfem2::opengl::DrawMeshTri3D_Edge(aXYZ, aTri);
      delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(aXYZ, aTri, aFlgElm, aColor);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
