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
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/srchuni_v3.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"
#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mat4.h"
#include "delfem2/thread.h"
#include "delfem2/sampling.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/tex.h"

namespace dfm2 = delfem2;

// ----------------------------------------

bool Intersection(
    dfm2::PointOnSurfaceMesh<double> &pos_mesh,
    const dfm2::CVec3d &src1,
    const dfm2::CVec3d &dir1,
    const std::vector<double> &vec_xyz,
    const std::vector<unsigned int> &vec_tri,
    const std::vector<dfm2::CNodeBVH2> &bvh_nodes,
    const std::vector<dfm2::CBV3_Sphere<double>> &bvh_volumes){
  std::vector<unsigned int> aIndElem;
  BVH_GetIndElem_Predicate(
      aIndElem,
      dfm2::CIsBV_IntersectLine<dfm2::CBV3_Sphere<double>, double>(src1.p, dir1.p),
      0, bvh_nodes, bvh_volumes);
  if (aIndElem.empty()) { return false; } // no bv hit the ray
  std::map<double, dfm2::PointOnSurfaceMesh<double>> mapDepthPES;
  IntersectionRay_MeshTri3DPart(
      mapDepthPES,
      src1, dir1,
      vec_tri, vec_xyz, aIndElem, 1.0e-10);
  if (mapDepthPES.empty()) { return false; }
  pos_mesh = mapDepthPES.begin()->second;
  return true;
}

double Radiance_Retroreflection(
    std::array<unsigned short, 3> &Xi,
    const dfm2::CVec3d &src1,
    const dfm2::CVec3d &dir1,
    const std::vector<double> &vec_xyz,
    const std::vector<unsigned int> &vec_tri,
    const std::vector<dfm2::CNodeBVH2> &bvh_nodes,
    const std::vector<dfm2::CBV3_Sphere<double>> &bvh_volumes,
    const dfm2::CVec3d &dir_view,
    int idepth){
  if( idepth >= 3 ){ return 0; }
  dfm2::PointOnSurfaceMesh<double> pos_mesh;
  bool is_hit = Intersection(
      pos_mesh,
      src1, dir1, vec_xyz, vec_tri, bvh_nodes, bvh_volumes);
  if( !is_hit ){ return 0; }
  unsigned int itri = pos_mesh.itri;
  assert(itri < vec_tri.size() / 3);
  double n[3], area;
  delfem2::UnitNormalAreaTri3(
      n, area,
      vec_xyz.data() + vec_tri[itri * 3 + 0] * 3,
      vec_xyz.data() + vec_tri[itri * 3 + 1] * 3,
      vec_xyz.data() + vec_tri[itri * 3 + 2] * 3);
  dfm2::CVec3d src2 = pos_mesh.PositionOnMeshTri3(vec_xyz, vec_tri);
  double d0 = -dir_view.dot(dfm2::CVec3d(n));
  dfm2::CVec3d src2a = src2 + dfm2::CVec3d(n) * 1.0e-3;
  double visiblity = 0.;
  if( idepth == 0 ) {
    visiblity = 1.0;
  }
  else if( d0 > 0 ){
    dfm2::PointOnSurfaceMesh<double> pos_mesh2;
    bool is_hit2 = Intersection(
        pos_mesh2,
        src2a, -dir_view, vec_xyz, vec_tri, bvh_nodes, bvh_volumes);
    if( !is_hit2 ){ visiblity = 1.; }
  }
  const dfm2::CVec3d dir2 = SampleHemisphereNormalCos(dfm2::CVec3d(n), delfem2::RandomVec2<double>(Xi));
  return d0 * visiblity + Radiance_Retroreflection(
      Xi,
      src2a,dir2,vec_xyz,vec_tri,bvh_nodes,bvh_volumes,dir_view,idepth+1);
}


double Radiance_AmbientOcclusion(
    std::array<unsigned short, 3> &Xi,
    const dfm2::CVec3d &src1,
    const dfm2::CVec3d &dir1,
    const std::vector<double> &vec_xyz,
    const std::vector<unsigned int> &vec_tri,
    const std::vector<dfm2::CNodeBVH2> &bvh_nodes,
    const std::vector<dfm2::CBV3_Sphere<double>> &bvh_volumes ) {
  dfm2::PointOnSurfaceMesh<double> pos_mesh;
  bool is_hit = Intersection(
      pos_mesh,
      src1, dir1, vec_xyz, vec_tri, bvh_nodes, bvh_volumes);
  if( !is_hit ){ return 0; }
  unsigned int itri = pos_mesh.itri;
  assert(itri < vec_tri.size() / 3);
  double n[3], area;
  delfem2::UnitNormalAreaTri3(
      n, area,
      vec_xyz.data() + vec_tri[itri * 3 + 0] * 3,
      vec_xyz.data() + vec_tri[itri * 3 + 1] * 3,
      vec_xyz.data() + vec_tri[itri * 3 + 2] * 3);
  dfm2::CVec3d src2 = pos_mesh.PositionOnMeshTri3(vec_xyz, vec_tri);
  dfm2::CVec3d src2a = src2 + dfm2::CVec3d(n) * 1.0e-3;
  const dfm2::CVec3d dir2 = SampleHemisphereNormalCos(dfm2::CVec3d(n), delfem2::RandomVec2<double>(Xi));
  dfm2::PointOnSurfaceMesh<double> pos_mesh2;
  bool is_hit2 = Intersection(
      pos_mesh2,
      src2a, dir2, vec_xyz, vec_tri, bvh_nodes, bvh_volumes);
  if( !is_hit2 ){ return 1; }
  return 0;
}

int main() {
  std::vector<double> vec_xyz; // 3d points
  std::vector<unsigned int> vec_tri;
  { // load input mesh
    delfem2::Read_Ply(
        vec_xyz, vec_tri,
        std::filesystem::path(PATH_INPUT_DIR) / "bunny_2k.ply");
    dfm2::Normalize_Points3(vec_xyz, 2.0);
  }
  std::vector<dfm2::CNodeBVH2> bvh_nodes;
  std::vector<dfm2::CBV3_Sphere<double>> bvh_volumes;
  delfem2::ConstructBVHTriangleMeshMortonCode(
      bvh_nodes, bvh_volumes,
      vec_xyz, vec_tri);
  const unsigned int nw = 256;
  const unsigned int nh = 256;
  // above: constant
  // -------------------------------
  // below: changing during execution
  std::vector<float> afRGB(nw * nh * 3, 0.f);
  unsigned int isample = 0;
  dfm2::CMat4d mMVPd_inv;
  unsigned int imode = 0;
  auto render = [&](int iw, int ih) {
    std::array<unsigned short, 3> Xi = {
        (unsigned short) (ih*ih),
        (unsigned short) (iw*iw),
        (unsigned short) (isample * isample) };
    const double ps[4] = {-1. + (2. / nw) * (iw + 0.5), -1. + (2. / nh) * (ih + 0.5), +1., 1.};
    const double pe[4] = {-1. + (2. / nw) * (iw + 0.5), -1. + (2. / nh) * (ih + 0.5), -1., 1.};
    std::array<double, 3> qs = mMVPd_inv.MultVec3_Homography(ps);
    std::array<double, 3> qe = mMVPd_inv.MultVec3_Homography(pe);
    const dfm2::CVec3d src1(qs.data());
    const dfm2::CVec3d dir1 = (dfm2::CVec3d(qe.data()) - src1).normalized();
    double rd0;
    if( imode == 0 ) {
      rd0 = Radiance_Retroreflection(
          Xi,
          src1, dir1, vec_xyz, vec_tri, bvh_nodes, bvh_volumes, dir1, 0);
    }
    else{
      rd0 = Radiance_AmbientOcclusion(
          Xi,
          src1, dir1, vec_xyz, vec_tri, bvh_nodes, bvh_volumes);
    }
    dfm2::CVec3d r_ave(rd0,rd0,rd0);
    {
      float *ptr = afRGB.data() + (ih * nw + iw) * 3;
      const auto isamplef = static_cast<float>(isample);
      ptr[0] = (isamplef * ptr[0] + r_ave[0]) / (isamplef + 1.f);
      ptr[1] = (isamplef * ptr[1] + r_ave[1]) / (isamplef + 1.f);
      ptr[2] = (isamplef * ptr[2] + r_ave[2]) / (isamplef + 1.f);
    }
  };

  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    tex.width = nw;
    tex.height = nh;
    tex.channels = 3;
    tex.pixel_color.resize(tex.width * tex.height * tex.channels);
  }

  dfm2::glfw::CViewer3 viewer(2.f);
  viewer.width = 400;
  viewer.height = 400;
  viewer.camerachange_callbacks.emplace_back(
      [&afRGB, &isample] {
        std::fill(afRGB.begin(), afRGB.end(), 0.0);
        isample = 0;
      }
  );
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  tex.InitGL();

  while (!glfwWindowShouldClose(viewer.window)) {
    {
      double t = glfwGetTime();
      unsigned int imode_old = imode;
      imode = static_cast<int>(t*0.5) % 2;
      if( imode != imode_old ){
        std::fill(afRGB.begin(), afRGB.end(), 0.0);
        isample = 0;
      }
    }
    {
      const dfm2::CMat4f mP = viewer.GetProjectionMatrix();
      const dfm2::CMat4f mMV = viewer.GetModelViewMatrix();
      const dfm2::CMat4d mMVP = (mP * mMV).cast<double>();
      mMVPd_inv = dfm2::Inverse_Mat4(mMVP.data());
    }
    dfm2::parallel_for(nw, nh, render);
    isample++;
    for (unsigned int ih = 0; ih < tex.height; ++ih) {
      for (unsigned int iw = 0; iw < tex.width; ++iw) {
        for (int ic = 0; ic < 3; ++ic) {
          float fc = afRGB[(ih * tex.width + iw) * 3 + ic];
          fc = (fc > 1.f) ? 1.f : fc;
          fc = (fc < 0.f) ? 0.f : fc;
          int ifc = int(fc * 255.f + .5f);
          tex.pixel_color[(ih * tex.width + iw) * 3 + ic] = ifc;
        }
      }
    }
    tex.InitGL();
    //
    viewer.DrawBegin_oldGL();
    ::glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    ::glClear(GL_COLOR_BUFFER_BIT);
    ::glDisable(GL_LIGHTING);
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
    tex.Draw_oldGL();
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


