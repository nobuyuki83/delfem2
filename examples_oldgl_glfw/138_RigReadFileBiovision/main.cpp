/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <vector>
#include <climits>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should be before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/rig_bvh.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/rigv3.h"

namespace dfm2 = delfem2;

// -------------------------------------------

int main() {
  std::vector<dfm2::CRigBone> aBone;
  std::vector<dfm2::CChannel_BioVisionHierarchy> aChannelRotTransBone;
  size_t nframe = 0;
  std::vector<double> aValRotTransBone;

  //std::string path_bvh = std::string(PATH_INPUT_DIR) + "/jump.bvh";
  std::string path_bvh = std::string(PATH_INPUT_DIR) + "/walk.bvh";
//  std::cout << "path:" << path_bvh << std::endl;
  {
    double frame_time;
    std::string header_bvh;
    Read_BioVisionHierarchy(
      aBone, aChannelRotTransBone, nframe, frame_time, aValRotTransBone, header_bvh,
      path_bvh);
  }
  UpdateBoneRotTrans(aBone);
  std::cout << "nBone:" << aBone.size() << "   aCh:" << aChannelRotTransBone.size() << std::endl;
  for (unsigned int ib = 0; ib < aBone.size(); ++ib) {
    std::cout << ib << " " << aBone[ib].name << std::endl;
  }

  // -------------------------

  delfem2::glfw::CViewer3 viewer(15.0);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  // -------------------------

  while (!glfwWindowShouldClose(viewer.window)) {
    {
      static int iframe = 0;
      const size_t nch = aChannelRotTransBone.size();
      SetPose_BioVisionHierarchy(
        aBone, aChannelRotTransBone,
        aValRotTransBone.data() + iframe * nch);
      iframe = (iframe + 1) % nframe;
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    delfem2::opengl::DrawAxis(10);
    /*
    dfm2::opengl::DrawBone_Line(
        aBone,
        -1, -1,
        0.1, 1.0);
        */
    dfm2::opengl::DrawBoneCurrent_Octahedron(
      aBone,
      UINT_MAX, UINT_MAX,
      0.1, 1.0);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}