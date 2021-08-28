/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/gizmo_geo3.h"
#include "delfem2/mshio.h"
#include "delfem2/points.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mat4.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// -------------------------------------------

int main(
    [[maybe_unused]] int argc,
    [[maybe_unused]] char* argv[])
{
  delfem2::glfw::CViewer3 viewer;
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Ply(
      std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
      aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<unsigned int> aTri_cage;
  std::vector<double> aXYZ_cage;
  {
    double pmin[3], pmax[3];
    delfem2::BoundingBox3_Points3(
        pmin, pmax,
        aXYZ.data(), aXYZ.size()/3);
    pmin[0] -= 0.1;
    pmin[1] -= 0.1;
    pmin[2] -= 0.1;
    pmax[0] += 0.1;
    pmax[1] += 0.1;
    pmax[2] += 0.1;
    std::vector<unsigned int> aQuad_cage;
    delfem2::MeshQuad3_CubeVox(
        aXYZ_cage, aQuad_cage,
        pmin, pmax);
    delfem2::convert2Tri_Quad(
        aTri_cage, aQuad_cage);
  }
  // --------------------
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  // --------------------
  while(true){
    viewer.DrawBegin_oldGL();
    ::glColor3d(0,0,0);
    delfem2::opengl::DrawMeshTri3D_Edge(
        aXYZ_cage.data(), aXYZ_cage.size()/3,
        aTri_cage.data(), aTri_cage.size()/3);
    delfem2::opengl::DrawMeshTri3D_Edge(
        aXYZ.data(), aXYZ.size()/3,
        aTri.data(), aTri.size()/3);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(
        aXYZ.data(),
        aTri.data(), aTri.size()/3);
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


