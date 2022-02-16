

#include <cstdio>
#include <iostream>
#include <vector>
#if defined(_MSC_VER)
#  include <windows.h>
#endif
//
#define GL_SILENCE_DEPRECATION
#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#  define GL_GLEXT_PROTOTYPES
#  define EGL_EGLEXT_PROTOTYPES
#else
#  include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

#include "delfem2/vec3.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_unindexed.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/srch_bruteforce.h"
#include "delfem2/srch_bvh.h"
#include "delfem2/srch_trimesh3_class.h"
#include "delfem2/opengl/new/drawer_mshunindex.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

// ---------------------------------------------------------------

class MyViewer : public delfem2::glfw::CViewer3
{
 public:
  void InitGL(){
    drawer.InitGL();
  }
  void SetMesh(
    const std::vector<double> &vtx_xyz_,
    const std::vector<unsigned int> &tri_vtx_){
    vtx_xyz = vtx_xyz_;
    tri_vtx = tri_vtx_;
    tri_flg.assign(tri_vtx.size() / 3, 0);
    {
      std::vector<double> tri_xyz, tri_nrm;
      dfm2::UnidexedVertexDataTriMesh(
        tri_xyz,
        vtx_xyz, tri_vtx);
      dfm2::UnindexedNormalTriMesh3(
        tri_nrm,
        vtx_xyz, tri_vtx);
      dfm2::UnindexedColorTriMesh3(
        tri_rgb,
        flg_rgb, tri_flg, tri_vtx);
      drawer.Initialize(tri_xyz, tri_nrm, tri_rgb, 3);
    }
    dfm2::ConstructBVHTriangleMeshMortonCode(
      bvhnodes, bvhnode_aabb,
      vtx_xyz, tri_vtx);
  }
  void draw(){
    ::glEnable(GL_DEPTH_TEST);
    drawer.Draw(
      this->GetProjectionMatrix().data(),
      this->GetModelViewMatrix().data());
  }
  void mouse_press(const float src[3], const float dir[3]) override {
    unsigned int itri = this->PickTri(src, dir);
    if (itri == UINT_MAX) { return; }
    tri_flg[itri] = 1;
    dfm2::UnindexedColorTriMesh3(
      tri_rgb,
      flg_rgb, tri_flg, tri_vtx);
    drawer.UpdateTriRgb(tri_rgb);
  }
 private:
  unsigned int PickTri(const float src[3], const float dir[3]) const {
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
    return mapDepthPES.begin()->second.itri;
  }
 public:
  dfm2::opengl::Drawer_MeshUnIndexed drawer;
  //
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  std::vector<int> tri_flg;
  std::vector<dfm2::CNodeBVH2> bvhnodes;
  std::vector<dfm2::CBV3_Sphere<double>> bvhnode_aabb;
  std::vector<double> tri_rgb;
  const double flg_rgb[2][3] = {
    {1, 1, 1},
    {1, 0.4, 0.4}};
};

int main(int, char **) {
  MyViewer viewer;
  //
  delfem2::glfw::InitGLNew();
  viewer.OpenWindow();
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  viewer.InitGL();
  {
    std::vector<double> vtx_xyz;
    std::vector<unsigned int> tri_vtx;
    std::string filePathName = std::string(PATH_SOURCE_DIR) + "/../../test_inputs/bunny_1k.obj";
    delfem2::Read_Obj3(
      filePathName,
      vtx_xyz, tri_vtx);
    delfem2::Normalize_Points3(
      vtx_xyz,
      1.);
    viewer.SetMesh(vtx_xyz, tri_vtx);
  }

  while (!glfwWindowShouldClose(viewer.window)) {
    int display_w, display_h;
    glfwGetFramebufferSize(viewer.window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    ::glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(1.1f, 4.0f);
    viewer.draw();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  return 0;
}
