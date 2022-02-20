/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <vector>
#include <set>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/pbd_geo3.h"
#include "delfem2/pbd_geo3dtri23.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/cad2.h"
#include "delfem2/cad2_mesher.h"
#include "delfem2/dtri.h"
#include "delfem2/mshmisc.h"
#include "delfem2/srch_trimesh3_class.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/tinygltf/io_gltf.h"
#include "tinygltf/tiny_gltf.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/caddtri_v3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// --------------------------------------------

std::vector<dfm2::CDynPntSur> aPo2D;
std::vector<dfm2::CDynTri> aETri;
std::vector<dfm2::CVec2d> aVec2;
std::vector<unsigned int> aLine;
std::vector<dfm2::CInfoNearest<double>> aInfoNearest;

dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> bvh;

const double dt = 0.01;
const double gravity[3] = {0.0, 0.0, 0.0};
const double contact_clearance = 0.0001;
const double rad_explore = 0.1;

// ------------------------------------

void StepTime(
    std::vector<double> &vtx_xyz_cloth,
    std::vector<double> &vtx_xyz_clothtmp,
    std::vector<double> &vtx_uvw_cloth,
    const std::vector<int> &vtx_bcflag,
    const std::vector<double> &vtx_xyz_body,
    const std::vector<unsigned int> &tri_vtx_body,
    const std::vector<double> &vtx_nrm_body) {
  dfm2::PBD_Pre3D(
      vtx_xyz_clothtmp,
      dt, gravity, vtx_xyz_cloth, vtx_uvw_cloth, vtx_bcflag);
  dfm2::PBD_TriStrain(
      vtx_xyz_clothtmp.data(),
      vtx_xyz_clothtmp.size() / 3, aETri, aVec2);
  dfm2::PBD_Bend(
      vtx_xyz_clothtmp.data(),
      vtx_xyz_clothtmp.size() / 3, aETri, aVec2, 1.0);
  dfm2::PBD_Seam(
      vtx_xyz_clothtmp.data(),
      vtx_xyz_clothtmp.size() / 3, aLine.data(), aLine.size() / 2);
  dfm2::Project_PointsIncludedInBVH_Outside_Cache(
      vtx_xyz_clothtmp.data(), aInfoNearest,
      vtx_xyz_clothtmp.size() / 3,
      contact_clearance, bvh,
      vtx_xyz_body.data(), vtx_xyz_body.size() / 3,
      tri_vtx_body.data(), tri_vtx_body.size() / 3,
      vtx_nrm_body.data(), rad_explore);
  dfm2::PBD_Post(
      vtx_xyz_cloth, vtx_uvw_cloth,
      dt, vtx_xyz_clothtmp, vtx_bcflag);
}

// ---------------------------------------

void myGlutDisplay(
    const std::vector<double> &vtx_xyz_cloth,
    const std::vector<double> &vtx_xyz_body,
    const std::vector<unsigned int> &tri_vtx_body) {
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  //  ::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1, 4.0);

  ::glPointSize(5);
  ::glLineWidth(1);
  {
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.8, 0.8, 0.8);
    /*
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    ::glEnable(GL_DEPTH_TEST);
     */
//    DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  }

  ::glDisable(GL_LIGHTING);
  ::glColor3d(0, 0, 0);
  delfem2::opengl::DrawMeshDynTri3D_Edge(vtx_xyz_cloth, aETri);

  ::glColor3d(1, 0, 0);
  delfem2::opengl::DrawMeshTri3D_Edge(
      vtx_xyz_body.data(), vtx_xyz_body.size() / 3,
      tri_vtx_body.data(), tri_vtx_body.size() / 3);
//  DrawSphere_Edge(rad0);
}

class CRigidTrans_2DTo3D {
 public:
  dfm2::CVec2d org2;
  dfm2::CVec3d org3;
  dfm2::CMat3d R;
};

int main() {
  delfem2::CCad2D cad;
  {
    const double xys0[8] = {-0.0, -0.0, +1.0, -0.0, +1.0, +1.0, -0.0, +1.0,};
    const double xys1[8] = {+2.0, -0.0, +3.0, -0.0, +3.0, +1.0, +2.0, +1.0,};
    cad.AddPolygon(std::vector<double>(xys0, xys0 + 8));
    cad.AddPolygon(std::vector<double>(xys1, xys1 + 8));
  }
  delfem2::CMesher_Cad2D mesher;
  mesher.edge_length = 0.04;
  delfem2::CMeshDynTri2D dmesh;
  mesher.Meshing(dmesh,
                 cad);
  dmesh.Check();
  aPo2D = dmesh.aEPo;
  aETri = dmesh.aETri;
  aVec2 = dmesh.aVec2;

  // -----------------------------
  std::vector<double> vtx_xyz_cloth; // deformed vertex positions
  std::vector<double> vtx_xyz_clothtmp;
  std::vector<double> vtx_uvw_cloth; // deformed vertex velocity
  std::vector<int> vtx_bcflag;  // boundary condition flag (0:free 1:fixed)
  const int np = aPo2D.size();
  vtx_uvw_cloth.resize(np * 3, 0.0);
  vtx_bcflag.resize(np, 0);
  vtx_xyz_cloth.resize(np * 3);
  {
    CRigidTrans_2DTo3D rt23;
    rt23.org2 = dfm2::CVec2d(2.5, 0.5);
    rt23.org3 = dfm2::CVec3d(0.0, 0.0, 0.5);
    rt23.R = dfm2::Mat3_RotMatFromAxisAngleVec(std::array<double, 3>{0.0, 3.1415, 0.0});
    std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1, 1), cad);
    for (int ip: aIP) {
      dfm2::CVec3d p0(aVec2[ip].x - rt23.org2.x, aVec2[ip].y - rt23.org2.y, 0.0);
      dfm2::CVec3d p1 = rt23.org3 + rt23.R * p0;
      vtx_xyz_cloth[ip * 3 + 0] = p1.x;
      vtx_xyz_cloth[ip * 3 + 1] = p1.y;
      vtx_xyz_cloth[ip * 3 + 2] = p1.z;
    }
    {
      CRigidTrans_2DTo3D rt23;
      rt23.org2 = dfm2::CVec2d(0.5, 0.5);
      rt23.org3 = dfm2::CVec3d(0.0, 0.0, -0.5);
      rt23.R.SetIdentity();
      std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1, 0), cad);
      for (int ip: aIP) {
        dfm2::CVec3d p0(aVec2[ip].x - rt23.org2.x, aVec2[ip].y - rt23.org2.y, 0.0);
        dfm2::CVec3d p1 = rt23.org3 + rt23.R * p0;
        vtx_xyz_cloth[ip * 3 + 0] = p1.x;
        vtx_xyz_cloth[ip * 3 + 1] = p1.y;
        vtx_xyz_cloth[ip * 3 + 2] = p1.z;
      }
    }
    aLine.clear();
    {
      std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(1, true, cad);
      std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(7, true, cad);
      const unsigned int npe = aIP0.size();
      assert(aIP1.size() == npe);
      for (unsigned int iip = 0; iip < npe; ++iip) {
        int ip0 = aIP0[iip];
        int ip1 = aIP1[npe - iip - 1];
        aLine.push_back(ip0);
        aLine.push_back(ip1);
      }
    }
    {
      std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(3, true, cad);
      std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(5, true, cad);
      const std::size_t npe = aIP0.size();
      assert(aIP1.size() == npe);
      for (unsigned int iip = 0; iip < npe; ++iip) {
        unsigned int ip0 = aIP0[iip];
        unsigned int ip1 = aIP1[npe - iip - 1];
        aLine.push_back(ip0);
        aLine.push_back(ip1);
      }
    }
  }
  vtx_xyz_clothtmp = vtx_xyz_cloth;

  std::vector<double> vtx_xyz_body;
  std::vector<unsigned int> tri_vtx_body;
  std::vector<double> vtx_nrm_body;
  { // make a unit sphere
    {
      tinygltf::Model model;
      tinygltf::TinyGLTF loader;
      std::string err;
      std::string warn;
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Monster.glb";
      std::string path_gltf = std::string(PATH_INPUT_DIR) + "/CesiumMan.glb";
      bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, path_gltf); // for binary glTF(.glb)
      if (!warn.empty()) { printf("Warn: %s\n", warn.c_str()); }
      if (!err.empty()) { printf("Err: %s\n", err.c_str()); }
      if (!ret) {
        printf("Failed to parse glTF\n");
        return -1;
      }
      dfm2::Print(model);
      //
      std::vector<double> vtx_xyz_bodyini;
      std::vector<double> skinning_sparse_weight;
      std::vector<unsigned int> skinning_sparse_index;
      std::vector<dfm2::CRigBone> bones;
      dfm2::GetMeshInfo(
          vtx_xyz_bodyini, tri_vtx_body,
          skinning_sparse_weight, skinning_sparse_index,
          model, 0, 0);
      dfm2::GetBone(bones, model, 0);
      {
        dfm2::Quat_Bryant(bones[0].quatRelativeRot, -90.0 / 180.0 * 3.1415, 0.0, 0.0);
        bones[0].transRelative[2] -= 1.0;
//        Quat_Bryant(aBone[0].rot, 0, 0, 0);
      }
      UpdateBoneRotTrans(bones);
      vtx_xyz_body = vtx_xyz_bodyini;
      dfm2::Skinning_LBS_LocalWeight(
          vtx_xyz_body.data(),
          vtx_xyz_bodyini.data(), vtx_xyz_bodyini.size() / 3,
          bones, skinning_sparse_weight.data(), skinning_sparse_index.data());
    }

    vtx_nrm_body.resize(vtx_xyz_body.size());
    delfem2::Normal_MeshTri3D(
        vtx_nrm_body.data(),
        vtx_xyz_body.data(), vtx_xyz_body.size() / 3,
        tri_vtx_body.data(), tri_vtx_body.size() / 3);
    bvh.Init(
        vtx_xyz_body.data(), vtx_xyz_body.size() / 3,
        tri_vtx_body.data(), tri_vtx_body.size() / 3,
        0.01);
  }

  delfem2::glfw::CViewer3 viewer;
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  // Enter main loop
  while (true) {
    StepTime(
        vtx_xyz_cloth, vtx_xyz_clothtmp, vtx_uvw_cloth, vtx_bcflag,
        vtx_xyz_body, tri_vtx_body, vtx_nrm_body);
    // ------------
    viewer.DrawBegin_oldGL();
    myGlutDisplay(
        vtx_xyz_cloth,
        vtx_xyz_body, tri_vtx_body);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
