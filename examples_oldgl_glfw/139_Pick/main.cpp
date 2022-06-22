/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should be before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/msh_affine_transformation.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/msh_center_of_gravity.h"
#include "delfem2/msh_boundingbox.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/srch_bvh.h"
#include "delfem2/srch_bruteforce.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"

namespace dfm2 = delfem2;

/*
class MyView
    : public dfm2::glfw::CViewer3 {
 public:

  MyView(
      double view_height,
      std::vector<unsigned int> &aFlagElem0,
      const std::vector<double> &aXYZ0, // 3d points
      const std::vector<unsigned int> &aTri0)
      : CViewer3(view_height), 
      tri_flg(aFlagElem0), vtx_xyz(aXYZ0), tri_vtx(aTri0) {

  }

  void mouse_press(const float src[3], const float dir[3]) override {
    unsigned int itri = this->PickTri(src, dir);
    if (itri == UINT_MAX) { return; }
    tri_flg[itri] = 1;
  }
 private:

 public:
  std::vector<unsigned int> &tri_flg;
  const std::vector<double> &vtx_xyz; // 3d points
  const std::vector<unsigned int> &tri_vtx;
};
 */

int main() {
  std::vector<double> vtx_xyz; // 3d points
  std::vector<unsigned int> tri_vtx;

  { // load input mesh
    delfem2::Read_Ply(
        vtx_xyz, tri_vtx,
        std::filesystem::path(PATH_SOURCE_DIR) / ".." / ".." / "test_inputs" / "arm_16k.ply");
    dfm2::Normalize_Points3(vtx_xyz, 2.0);
    std::cout << "point_size: " << vtx_xyz.size() / 3 << std::endl;
    std::cout << "triangle_size: " << tri_vtx.size() / 3 << std::endl;
  }

  std::vector<unsigned int> aFlagElem(tri_vtx.size() / 3, 0);
  std::vector<dfm2::CNodeBVH2> bvhnodes;
  std::vector<dfm2::CBV3_Sphere<double>> bvhnode_aabb;
  std::vector<unsigned int> tri_adjtri;
  {
    { // make BVH
      std::vector<double> aCent;
      dfm2::CentsMaxRad_MeshTri3(
          aCent,
          vtx_xyz, tri_vtx);
      double min_xyz[3], max_xyz[3];
      delfem2::BoundingBox3_Points3(
          min_xyz, max_xyz,
          aCent.data(),
          aCent.size() / 3);
      std::vector<unsigned int> aSortedId;
      std::vector<std::uint32_t> aSortedMc;
      dfm2::SortedMortenCode_Points3(
          aSortedId, aSortedMc,
          aCent, min_xyz, max_xyz);
      dfm2::BVHTopology_Morton(
          bvhnodes,
          aSortedId, aSortedMc);
      dfm2::CLeafVolumeMaker_Mesh<dfm2::CBV3_Sphere<double>, double> lvm(
          1.0e-10,
          vtx_xyz.data(), vtx_xyz.size() / 3,
          tri_vtx.data(), tri_vtx.size() / 3, 3);
      dfm2::BVH_BuildBVHGeometry(
          bvhnode_aabb,
          0, bvhnodes,
          lvm);
      { // assertion
        dfm2::Check_MortonCode_Sort(
            aSortedId, aSortedMc, aCent, min_xyz, max_xyz);
        dfm2::Check_MortonCode_RangeSplit(aSortedMc);
        dfm2::Check_BVH(
            bvhnodes,
            aCent.size() / 3);
      }
    }
    { // make triangle surrounding graph
      std::vector<unsigned int> elsup_ind, elsup;
      dfm2::JArray_ElSuP_MeshElem(
          elsup_ind, elsup,
          tri_vtx.data(), tri_vtx.size() / 3, 3,
          vtx_xyz.size() / 3);
      ElSuEl_MeshElem(
          tri_adjtri,
          tri_vtx.data(), tri_vtx.size() / 3,
          delfem2::MESHELEM_TRI,
          vtx_xyz.size() / 3);
    }
  }

  delfem2::glfw::CViewer3 viewer(1.5);
  viewer.mousepress_callbacks.emplace_back(
      [&aFlagElem,&bvhnodes,&bvhnode_aabb,&tri_vtx,&vtx_xyz](const float src[3], const float dir[3]){
        std::vector<unsigned int> aIndElem;
        dfm2::BVH_GetIndElem_Predicate(
            aIndElem,
            delfem2::CIsBV_IntersectLine<dfm2::CBV3_Sphere<double>, float>(src, dir),
            0, bvhnodes, bvhnode_aabb);
        std::map<double, dfm2::PointOnSurfaceMesh<double>> mapDepthPES;
        dfm2::IntersectionRay_MeshTri3DPart(
            mapDepthPES,
            dfm2::CVec3d(src), dfm2::CVec3d(dir),
            tri_vtx, vtx_xyz, aIndElem,
            1.0e-3);
        if (mapDepthPES.empty()) { return UINT_MAX; }
        unsigned int itri = mapDepthPES.begin()->second.itri;
        aFlagElem[itri] = 1;
      }
  );

  // ----------------------
  // opengl start from here
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  const std::vector<std::pair<int, delfem2::CColor> > aColor = {
      {2, delfem2::CColor::White()},
      {2, delfem2::CColor::Red()}};

  while (!glfwWindowShouldClose(viewer.window)) {
    // ----------
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshTri3D_Edge(vtx_xyz, tri_vtx);
    delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(
        vtx_xyz, tri_vtx, aFlagElem, aColor);
    // ----------
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


