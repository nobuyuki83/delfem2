/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <set>
#include <vector>
#include <string>
#include <cassert>
#include <random>
#include <cstdlib>
#include <set>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshuni.h"
#include "delfem2/vec3.h"
#include "delfem2/point_on_surface_mesh.h"
#include "delfem2/geoproximity3_v3.h"

class RandomSamplingOnSurfaceOnMeshTri3 {
 public:
  RandomSamplingOnSurfaceOnMeshTri3(
    std::vector<double> &vtx_xyz_,
    std::vector<unsigned int> &tri_vtx_)
    : vtx_xyz(vtx_xyz_), tri_vtx(tri_vtx_) {
    unsigned int ntri = tri_vtx.size() / 3;
    cumulative_area_sum.reserve(ntri + 1);
    cumulative_area_sum.push_back(0.);
    for (unsigned int itri = 0; itri < ntri; ++itri) {
      unsigned int i0 = tri_vtx[itri * 3 + 0];
      unsigned int i1 = tri_vtx[itri * 3 + 1];
      unsigned int i2 = tri_vtx[itri * 3 + 2];
      double a0 = delfem2::Area_Tri3(
        vtx_xyz.data() + i0 * 3,
        vtx_xyz.data() + i1 * 3,
        vtx_xyz.data() + i2 * 3);
      double t0 = cumulative_area_sum[cumulative_area_sum.size() - 1];
      cumulative_area_sum.push_back(a0 + t0);
    }
    rndeng.seed(std::random_device{}());
    dist_01 = std::uniform_real_distribution<double>(0, 1);
  }
  std::tuple<unsigned int, double, double> Sample() {
    double val01 = dist_01(rndeng);
    unsigned int ntri = tri_vtx.size() / 3;
    assert(cumulative_area_sum.size() == ntri + 1);
    double a0 = val01 * cumulative_area_sum[ntri];
    unsigned int itri_l = 0, itri_u = ntri;
    for (;;) {
      assert(cumulative_area_sum[itri_l] < a0);
      assert(a0 < cumulative_area_sum[itri_u]);
      unsigned int itri_h = (itri_u + itri_l) / 2;
      if (itri_u - itri_l == 1) { break; }
      if (cumulative_area_sum[itri_h] < a0) { itri_l = itri_h; }
      else { itri_u = itri_h; }
    }
    assert(cumulative_area_sum[itri_l] < a0);
    assert(a0 < cumulative_area_sum[itri_l + 1]);
    double r0 = (a0 - cumulative_area_sum[itri_l]) / (cumulative_area_sum[itri_l + 1] - cumulative_area_sum[itri_l]);
    double r1 = dist_01(rndeng);
    if (r0 + r1 > 1) {
      double r0a = r0, r1a = r1;
      r0 = 1 - r1a;
      r1 = 1 - r0a;
    }
    return {itri_l, r0, r1};
  }
 public:
  const std::vector<double> &vtx_xyz;
  const std::vector<unsigned int> &tri_vtx;
  std::vector<double> cumulative_area_sum;
  std::mt19937 rndeng;
  std::uniform_real_distribution<double> dist_01;
};

std::vector<unsigned int> NeighbourElementIndexes(
  const std::array<double, 3> pos,
  double rad,
  unsigned int itri0,
  const std::vector<double> &vtx_xyz,
  const std::vector<unsigned int> &tri_vtx,
  const std::vector<unsigned int> &tri_adjtri) {
  std::set<unsigned int> buff;
  std::stack<unsigned int> st;
  st.push(itri0);
  while (!st.empty()) {
    unsigned int iel0 = st.top();
    st.pop();
    if (buff.find(iel0) != buff.end()) { continue; } // already studied
    double dist_min = -1;
    {
      double pn[3], r0, r1;
      delfem2::GetNearest_TrianglePoint3D(
        pn,r0,r1,
        pos.data(),
        vtx_xyz.data() + tri_vtx[iel0 * 3 + 0] * 3,
        vtx_xyz.data() + tri_vtx[iel0 * 3 + 1] * 3,
        vtx_xyz.data() + tri_vtx[iel0 * 3 + 2] * 3 );
      dist_min = delfem2::Distance3(pn,pos.data());
    }
    if (dist_min > rad) { continue; }
    buff.insert(iel0);
    for (unsigned int ie = 0; ie < 3; ++ie) {
      unsigned int iel1 = tri_adjtri[iel0 * 3 + ie];
      if (iel1 == UINT_MAX) { continue; }
      st.push(iel1);
    }
  }
  return {buff.begin(), buff.end()};
}

// ===========

void Draw(
  const std::vector<double> &vtx_xyz,
  const std::vector<unsigned int> &tri_vtx,
  const std::vector<std::tuple<unsigned int, double, double> > &samples) {
  ::glEnable(GL_LIGHTING);
  delfem2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz, tri_vtx);
  ::glDisable(GL_LIGHTING);
  ::glPointSize(5);
  ::glBegin(GL_POINTS);
  for (const auto &s: samples) {
    delfem2::PointOnSurfaceMesh<double> posm{
      std::get<0>(s), std::get<1>(s), std::get<2>(s)};
    auto pos = posm.PositionOnMeshTri3(vtx_xyz, tri_vtx);
    ::glVertex3d(pos[0], pos[1], pos[2]);
  }
  ::glEnd();
}

int main() {
  namespace dfm2 = delfem2;
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  dfm2::Read_Ply(
    vtx_xyz, tri_vtx,
    std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.ply");
  dfm2::Normalize_Points3(vtx_xyz);
  std::vector<unsigned int> tri_adjtri;
  dfm2::ElSuEl_MeshElem(
    tri_adjtri,
    tri_vtx.data(), tri_vtx.size() / 3,
    dfm2::MESHELEM_TRI, vtx_xyz.size() / 3);
  RandomSamplingOnSurfaceOnMeshTri3 mapper(vtx_xyz, tri_vtx);
  delfem2::glfw::CViewer3 viewer(1.0);
  // ----------
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  while (!glfwWindowShouldClose(viewer.window)) {
    std::vector<std::tuple<unsigned int, double, double> > samples;
    for (unsigned int ismpl = 0; ismpl < 10000; ++ismpl) {
      samples.push_back(mapper.Sample());
      viewer.DrawBegin_oldGL();
      ::glColor3d(1, 0, 0);
      Draw(vtx_xyz, tri_vtx, samples);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if( glfwWindowShouldClose(viewer.window) ){ break; }
    }
    if( glfwWindowShouldClose(viewer.window) ){ break; }
    //
    samples.clear();
    std::multimap<unsigned int, unsigned int> el_smpl;
    double rad = 0.05;
    for (unsigned int itr = 0; itr < 10000; ++itr) {
      const auto smpli = mapper.Sample();
      bool is_near = false;
      // check neighbour
      auto posi = delfem2::PointOnSurfaceMesh<double>{
        std::get<0>(smpli), std::get<1>(smpli), std::get<2>(smpli)
      }.PositionOnMeshTri3(vtx_xyz, tri_vtx);
      std::vector<unsigned int> aIE = NeighbourElementIndexes(
        posi, rad,
        std::get<0>(smpli), vtx_xyz, tri_vtx, tri_adjtri);
      for (auto ie: aIE) {
        auto[il, iu] = el_smpl.equal_range(ie);
        for (auto it = il; it != iu; ++it) {
          unsigned int jsmpl = it->second;
          const auto smplj = samples[jsmpl];
          const auto posj = delfem2::PointOnSurfaceMesh<double>{
            std::get<0>(smplj), std::get<1>(smplj), std::get<2>(smplj)
          }.PositionOnMeshTri3(vtx_xyz, tri_vtx);
          const double dist = dfm2::Distance3(posi.data(), posj.data());
          if (dist < rad) {
            is_near = true;
            break;
          }
        }
        if (is_near) { break; }
      }
      if (!is_near) {
        el_smpl.insert(std::make_pair(std::get<0>(smpli), samples.size()));
        samples.push_back(smpli);
      }
      //
      viewer.DrawBegin_oldGL();
      ::glColor3d(0, 0, 1);
      Draw(vtx_xyz, tri_vtx, samples);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if( glfwWindowShouldClose(viewer.window) ){ break; }
    }
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
