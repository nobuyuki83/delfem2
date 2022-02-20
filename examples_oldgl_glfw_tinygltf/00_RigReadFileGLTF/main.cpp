/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/tinygltf/io_gltf.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/rig_geo3.h"

namespace dfm2 = delfem2;

// ---------------------------------

int main() {
  std::vector<unsigned int> tri_vtx;
  std::vector<double> skinning_sparse_weight;
  std::vector<unsigned int> skinning_sparse_index;
  std::vector<double> vtx_xyz;
  std::vector<dfm2::CRigBone> bones;
  {
    //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
    //    std::string path_glb = std::string(PATH_INPUT_DIR)+"/Monster.glb";
    //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
    //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
    std::string path_glb = std::string(PATH_INPUT_DIR) + "/CesiumMan.glb";
    dfm2::CGLTF gltf;
    gltf.Read(path_glb);
    // gltf.Print();
    gltf.GetMeshInfo(
        vtx_xyz, tri_vtx, skinning_sparse_weight, skinning_sparse_index,
        0, 0);
    gltf.GetBone(bones, 0);
    // dfm2::SetCurrentBoneRotationAsDefault(bones);
  }
  const std::vector<double> vtx_xyz_ini = vtx_xyz;

  const dfm2::CQuatd quat_joint5_ini(bones[5].quatRelativeRot);
  {
    dfm2::CQuatd q0 = dfm2::Quat_Bryant(-0.5 * M_PI, -0.5 * M_PI, 0.0);
    q0 = q0 * dfm2::CQuatd(bones[0].quatRelativeRot);
    q0.CopyTo(bones[0].quatRelativeRot);
    bones[0].transRelative[0] = 1.0;
  }

  // --------------
  // opengl starts here
  delfem2::glfw::CViewer3 viewer(2);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window)) {
    iframe++;
    {
      const dfm2::CQuatd q0 = dfm2::Quat_Bryant(0.0, sin(iframe * 0.1), 0.0);
      (q0 * quat_joint5_ini).CopyTo(bones[5].quatRelativeRot);
      UpdateBoneRotTrans(bones);
      dfm2::Skinning_LBS_LocalWeight(
          vtx_xyz.data(),
          vtx_xyz_ini.data(), vtx_xyz_ini.size() / 3,
          bones, skinning_sparse_weight.data(), skinning_sparse_index.data());
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawAxis(1);
    // draw reference
    glColor3d(1,0,0);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vtx_xyz_ini.data(), vtx_xyz_ini.size() / 3,
        tri_vtx.data(), tri_vtx.size() / 3);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, std::array<float,3>{1,0,0}.data());
    delfem2::opengl::DrawBoneReference_Octahedron(
        bones,
        -1, -1,
        0.01, 1.0);
    // draw current
    glColor3d(0,0,1);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vtx_xyz.data(), vtx_xyz.size() / 3,
        tri_vtx.data(), tri_vtx.size() / 3);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, std::array<float,3>{0,0,1}.data());
    delfem2::opengl::DrawBoneCurrent_Octahedron(
        bones,
        -1, -1,
        0.01, 1.0);
    /*
    ::glDisable(GL_DEPTH_TEST);
    delfem2::opengl::DrawBone_Line(
        bones,
        -1, -1,
        0.01, 1.0);
    ::glEnable(GL_DEPTH_TEST);
     */
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


