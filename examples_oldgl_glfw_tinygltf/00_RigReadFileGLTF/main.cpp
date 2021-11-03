/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/tinygltf/io_gltf.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/rig_geo3.h"
#include <vector>
#include <set>
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ---------------------------------

int main()
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<double> aSkinSparseW;
  std::vector<unsigned int> aSkinSparseI;
  std::vector<double> aXYZ;
  std::vector<dfm2::CRigBone> aBone;
  {
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
//    std::string path_glb = std::string(PATH_INPUT_DIR)+"/Monster.glb";
    
    //      std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
    //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
    std::string path_glb = std::string(PATH_INPUT_DIR)+"/CesiumMan.glb";
    dfm2::CGLTF gltf;
    gltf.Read(path_glb);
    gltf.Print();
    gltf.GetMeshInfo(
        aXYZ0, aTri, aSkinSparseW, aSkinSparseI,
        0, 0);
    gltf.GetBone(aBone, 0);
    dfm2::SetCurrentBoneRotationAsDefault(aBone);
  }
  aXYZ = aXYZ0;

  // --------------
  // opengl starts here
  delfem2::glfw::CViewer3 viewer(2);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  int iframe = 0;
  while(!glfwWindowShouldClose(viewer.window)){
    iframe++;
    {
      dfm2::Quat_Bryant(aBone[0].quatRelativeRot,0.0,-0.5*M_PI,-0.5*M_PI);
      double dq[4]; dfm2::Quat_Bryant(dq,sin(iframe*0.1),0.0, 0.0);
      dfm2::Copy_Quat(aBone[5].quatRelativeRot, dq);
      UpdateBoneRotTrans(aBone);
      dfm2::Skinning_LBS_LocalWeight(
          aXYZ.data(),
          aXYZ0.data(), aXYZ0.size()/3,
          aBone, aSkinSparseW.data(), aSkinSparseI.data());
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ.data(), aTri.data(), aTri.size()/3);
    delfem2::opengl::DrawAxis(1);
    ::glDisable(GL_DEPTH_TEST);
    delfem2::opengl::DrawBone_Line(
        aBone,
        -1, -1,
        0.01, 1.0);
    ::glEnable(GL_DEPTH_TEST);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


