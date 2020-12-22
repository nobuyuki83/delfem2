/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief edit rigged shape
 * @details rotation for each bone and translation for root bone
 */

#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/opengl/glfw/ViewRig.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/tex.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/mat4.h"
#include <cstdlib>
#include <fstream>
#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace dfm2 = delfem2;

// ---------------------------

int main()
{
  dfm2::opengl::CViewerGlfw_RiggedMesh viewer;
  {
    std::vector<double> aXYZ0, aXYZ1;
    std::vector<unsigned int> aTri;
    std::vector<double> aWeightRigSparse;
    std::vector<unsigned int> aIdBoneRigSparse;
    std::vector<dfm2::CRigBone> aBone;
    //
    std::vector<double> aW;
    std::vector<unsigned int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl_Bone(
        aXYZ0,
        aW,
        aTri,
        aIndBoneParent,
        aJntRgrs,
        std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    {
      std::vector<double> aJntPos0;
      dfm2::Points3_WeighttranspPosition(
          aJntPos0,
          aJntRgrs, aXYZ0);
      dfm2::InitBones_JointPosition(
          aBone,
          aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
    }
    dfm2::SparsifyMatrixRow(
        aWeightRigSparse, aIdBoneRigSparse,
        aW.data(), aXYZ0.size()/3, aBone.size(),
        1.0e-5);
    viewer.SetRiggedMesh(aXYZ0,aTri,aWeightRigSparse,aIdBoneRigSparse,aBone);
  }
  // -----------
  dfm2::opengl::CTexRGB_Rect2D tex;
  int width, height;
  {
    int channels;
    unsigned char *img = stbi_load((std::string(PATH_INPUT_DIR)+"/uglysweater.jpg").c_str(),
                                   &width, &height, &channels, 0);
    tex.Initialize(width, height, img, "rgb");
    double scale = 0.0020;
    delete[] img;
    tex.max_x = -scale*width*0.5;
    tex.min_x = +scale*width*0.5;
    tex.max_y = -scale*height*0.5;
    tex.min_y = +scale*height*0.5;
    tex.z = -0.5;
    std::cout << width << " " << height << std::endl;
  }
  // -------------------
  viewer.Init_oldGL();
  tex.InitGL();
  dfm2::opengl::setSomeLighting();
  while ( !glfwWindowShouldClose(viewer.window) )
  {
    viewer.DrawBegin_oldGL();
    tex.Draw_oldGL();
    viewer.Draw();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
