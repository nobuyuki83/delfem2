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
#include "delfem2/rig_bvh.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/rig_geo3.h"

namespace dfm2 = delfem2;

// ---------------------------------

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  std::vector<double> skinning_sparse_weight;
  std::vector<unsigned int> skinning_sparse_index;
  std::vector<dfm2::CRigBone> bones_gltf;
  {
    std::string path_glb = std::string(PATH_SRC_DIR) + "/../../test_inputs/rig_transfer_from_bvh.glb";
    dfm2::CGLTF gltf;
    gltf.Read(path_glb);
    //gltf.Print();
    gltf.GetMeshInfo(
        vtx_xyz, tri_vtx, skinning_sparse_weight, skinning_sparse_index,
        0, 0);
    gltf.GetBone(bones_gltf, 0);
    // UpdateBoneRotTrans(bones_gltf);
    //dfm2::SetCurrentBoneRotationAsDefault(bones_gltf);
  }
  const std::vector<double> vtx_xyz_ini = vtx_xyz;
  // --------
  std::vector<dfm2::CRigBone> bones_bvh;
  std::vector<dfm2::CChannel_BioVisionHierarchy> aChannelRotTransBone;
  size_t nframe = 0;
  std::vector<double> aValRotTransBone;
  std::string path_bvh = std::string(PATH_SRC_DIR) + "/../../test_inputs/walk.bvh";
  {
    double frame_time;
    std::string header_bvh;
    Read_BioVisionHierarchy(
        bones_bvh, aChannelRotTransBone, nframe, frame_time, aValRotTransBone, header_bvh,
        path_bvh);
    UpdateBoneRotTrans(bones_bvh);
  }
  // std::cout << bones_gltf.size() << " " << bones_bvh.size() << std::endl;
  std::vector<unsigned int> mapGltf2Bvh(bones_gltf.size(), UINT_MAX);
  for (unsigned int igltf = 0; igltf < bones_gltf.size(); ++igltf) {
    // std::cout << igltf << " " << bones_gltf[igltf].name << std::endl;
    unsigned int jbvh = 0;
    for (jbvh = 0; jbvh < bones_bvh.size(); ++jbvh) {
      if (bones_gltf[igltf].name != bones_bvh[jbvh].name) { continue; }
      mapGltf2Bvh[igltf] = jbvh;
      break;
    }
    //delfem2::Print_Mat4(bones_gltf[igltf].affmat3Global);
    //delfem2::Print_Mat4(bones_gltf[igltf].invBindMat);
    //delfem2::CMat4d m = (delfem2::CMat4d(bones_gltf[igltf].affmat3Global)).Inverse();
    //delfem2::Print_Mat4(m);
    //delfem2::Print_Mat4(bones_bvh[jbvh].affmat3Global);
    //delfem2::Print_Mat4(bones_bvh[jbvh].invBindMat);
  }

  for (unsigned int igltf1 = 0; igltf1 < bones_gltf.size(); ++igltf1) {
    unsigned int jbvh1 = mapGltf2Bvh[igltf1];
    assert(jbvh1 < bones_bvh.size());
    dfm2::CVec3d p1 = dfm2::CMat4d(bones_bvh[jbvh1].invBindMat).Inverse().GetTranslationComponent();
    dfm2::CVec3d q1 = dfm2::CMat4d(bones_gltf[igltf1].invBindMat).Inverse().GetTranslationComponent();
    {
      const dfm2::CMat4d Tq1 = dfm2::CMat4d::Translation(-q1);
      const dfm2::CMat4d Tp1 = dfm2::CMat4d::Translation(p1);
      auto hoge = dfm2::CMat4d(bones_bvh[jbvh1].invBindMat) * Tp1 * Tq1;
      hoge.CopyTo(bones_gltf[igltf1].invBindMat);
    }
    {
      unsigned int jbvh0 = bones_bvh[jbvh1].ibone_parent;
      if (jbvh0 >= bones_bvh.size()) { continue; }
      unsigned int igltf0 = bones_gltf[igltf1].ibone_parent;
      if (igltf0 >= bones_gltf.size()) { continue; }
      dfm2::CVec3d p0 = dfm2::CMat4d(bones_bvh[jbvh0].invBindMat).Inverse().GetTranslationComponent();
      dfm2::CVec3d q0 = dfm2::CMat4d(bones_gltf[igltf0].invBindMat).Inverse().GetTranslationComponent();
      assert(bones_bvh[jbvh0].name == bones_gltf[igltf0].name);
      const dfm2::CMat4d Tq0 = dfm2::CMat4d::Translation(-q0);
      const dfm2::CMat4d Tp0 = dfm2::CMat4d::Translation(p0);
      dfm2::CMat4d R = dfm2::CMat4d::Identity();
      if ((p1 - p0).norm() > 1.0e-10) {
        R = dfm2::CMat4d::Mat3(Mat3_MinimumRotation(q1 - q0, p1 - p0).data());
      }
      auto hoge = dfm2::CMat4d(bones_bvh[jbvh0].invBindMat) * Tp0 * R * Tq0;
      hoge.CopyTo(bones_gltf[igltf0].invBindMat);
    }
  }

  // --------------
  // opengl starts here
  delfem2::glfw::CViewer3 viewer(2);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  unsigned int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window)) {
    iframe = (iframe + 1) % nframe;
    {
      const size_t nch = aChannelRotTransBone.size();
      SetPose_BioVisionHierarchy(
          bones_bvh, aChannelRotTransBone,
          aValRotTransBone.data() + iframe * nch);
      for (unsigned int igltf = 0; igltf < bones_gltf.size(); ++igltf) {
        unsigned int jbvh = mapGltf2Bvh[igltf];
        assert(jbvh<bones_bvh.size());
        delfem2::Copy_Mat4(
            bones_gltf[igltf].affmat3Global,
            bones_bvh[jbvh].affmat3Global);
        dfm2::Skinning_LBS_LocalWeight(
            vtx_xyz.data(),
            vtx_xyz_ini.data(), vtx_xyz_ini.size() / 3,
            bones_gltf, skinning_sparse_weight.data(), skinning_sparse_index.data());
      }
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    delfem2::opengl::DrawAxis(1);
    /*
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(
        vtx_xyz.data(), tri_vtx.data(), tri_vtx.size() / 3);
        */
    //
    //::glDisable(GL_DEPTH_TEST);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, std::array<float,4>{0,0,1}.data());
    delfem2::opengl::DrawBoneCurrent_Octahedron(
        bones_gltf,
        -1, -1,
        0.1, 1.0);
    glColor3fv(std::array<float,4>{0,0,1}.data());
    delfem2::opengl::DrawMeshTri3D_Edge(
        vtx_xyz.data(), vtx_xyz.size() /3,
        tri_vtx.data(), tri_vtx.size() / 3);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, std::array<float,4>{0.5,0.5,1}.data());
    delfem2::opengl::DrawBoneReference_Octahedron(
        bones_gltf,
        -1, -1,
        0.1, 1.0);
    delfem2::opengl::DrawBoneReference_Line(
        bones_gltf, 0.2);
    glColor3fv(std::array<float,4>{0.5,0.5,1}.data());
    delfem2::opengl::DrawMeshTri3D_Edge(
        vtx_xyz_ini.data(), vtx_xyz_ini.size() /3,
        tri_vtx.data(), tri_vtx.size() / 3);
    // ---------------------
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, std::array<float,4>{1,0,0}.data());
    delfem2::opengl::DrawBoneCurrent_Octahedron(
        bones_bvh,
        -1, -1,
        0.01, 1.0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, std::array<float,4>{1,0.5,0.5}.data());
    delfem2::opengl::DrawBoneReference_Octahedron(
        bones_bvh,
        -1, -1,
        0.01, 1.0);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


