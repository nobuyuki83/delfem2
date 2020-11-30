/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/tinygltf/io_gltf.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/rigv3_glold.h"
#include "delfem2/rig_geo3.h"
#include <vector>
#include <set>
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ---------------------------------

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<double> aRigWeight;
  std::vector<unsigned int> aRigJoint;
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
        aXYZ0, aTri, aRigWeight, aRigJoint,
        0,0);
    gltf.GetBone(aBone, 0);
  }
  aXYZ = aXYZ0;
  {
    double dq[4]; dfm2::Quat_Bryant(dq,0.5,0.,0.);
    double q0[4]; dfm2::Copy_Quat(q0, aBone[5].quatRelativeRot);
    dfm2::QuatQuat(aBone[5].quatRelativeRot,q0,dq);
    UpdateBoneRotTrans(aBone);
    dfm2::Skinning_LBS_LocalWeight(
        aXYZ.data(),
        aXYZ0.data(), aXYZ0.size()/3,
        aBone, aRigWeight.data(), aRigJoint.data());
  }
  // --------------
  // opengl starts here
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  while(!glfwWindowShouldClose(viewer.window)){
    // --------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ.data(), aTri.data(), aTri.size()/3);
    delfem2::opengl::DrawAxis(1);
    ::glDisable(GL_DEPTH_TEST);
    delfem2::opengl::DrawBone(
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


