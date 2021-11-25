/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <set>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/sampler_trimesh.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshuni.h"
#include "delfem2/point_on_surface_mesh.h"
#include "delfem2/geoproximity3_v3.h"
#include "delfem2/srchuni_v3.h"

// ===========

void Draw(
  const std::vector<double> &vtx_xyz,
  const std::vector<unsigned int> &tri_vtx,
  const std::vector<std::tuple<unsigned int, double, double> > &samples) {
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
  dfm2::RandomSamplingOnMeshTri3 mapper(vtx_xyz, tri_vtx);
  delfem2::glfw::CViewer3 viewer(1.0);
  // ----------
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  while (!glfwWindowShouldClose(viewer.window)) {
    {  // monte-carlo sampling on a triangle mesh
      std::vector<std::tuple<unsigned int, double, double> > samples;
      dfm2::RandomSamplingOnMeshTri3 mapper(vtx_xyz, tri_vtx);
      for (unsigned int ismpl = 0; ismpl < 3000; ++ismpl) {
        samples.push_back(mapper.Sample());
        viewer.DrawBegin_oldGL();
        ::glEnable(GL_LIGHTING);
        delfem2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz, tri_vtx);
        ::glColor3d(1, 0, 0);
        Draw(vtx_xyz, tri_vtx, samples);
        glfwSwapBuffers(viewer.window);
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
      if (glfwWindowShouldClose(viewer.window)) { break; }
    }
    {  // poisson sampling on triangle mesh
      std::vector<std::tuple<unsigned int, double, double> > samples;
      dfm2::RandomSamplingOnMeshTri3 mapper(vtx_xyz, tri_vtx);
      std::multimap<unsigned int, unsigned int> el_smpl;
      double rad = 0.05;
      for (unsigned int itr = 0; itr < 10000; ++itr) {
        const auto smpli = mapper.Sample();
        bool is_near = false;
        // check neighbour
        const auto posi = delfem2::PointOnSurfaceMesh<double>{
          std::get<0>(smpli), std::get<1>(smpli), std::get<2>(smpli)
        }.PositionOnMeshTri3(vtx_xyz, tri_vtx);
        std::vector<unsigned int> aIE = dfm2::IndexesOfConnectedTriangleInSphere(
          posi, rad,
          std::get<0>(smpli), vtx_xyz, tri_vtx, tri_adjtri);
        for (auto ie: aIE) {
          const auto[il, iu] = el_smpl.equal_range(ie);
          for (auto it = il; it != iu; ++it) {
            const unsigned int jsmpl = it->second;
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
        ::glEnable(GL_LIGHTING);
        delfem2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz, tri_vtx);
        ::glColor3d(0, 0, 1);
        Draw(vtx_xyz, tri_vtx, samples);
        glfwSwapBuffers(viewer.window);
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
    }
    {  // monte-carlo sampling on a triangle mesh
      std::vector<int> tri_flg(tri_vtx.size()/3,0);
      for(unsigned int it=0;it<tri_vtx.size()/3;++it){
        const double* p0 = vtx_xyz.data() + tri_vtx[it * 3 + 0] * 3;
        const double* p1 = vtx_xyz.data() + tri_vtx[it * 3 + 1] * 3;
        const double* p2 = vtx_xyz.data() + tri_vtx[it * 3 + 2] * 3;
        double y0 = (p0[1] + p1[1] + p2[1])/3.0;
        if( y0 > 0 ){ tri_flg[it] = 1; }
      }
      std::vector<std::tuple<unsigned int, double, double> > samples;
      dfm2::RandomSamplingOnMeshTri3Selective mapper(
        vtx_xyz, tri_vtx,
        [&tri_flg](unsigned int itri){ return tri_flg[itri] == 0; });
      for (unsigned int ismpl = 0; ismpl < 3000; ++ismpl) {
        samples.push_back(mapper.Sample());
        viewer.DrawBegin_oldGL();
        ::glEnable(GL_LIGHTING);
        delfem2::opengl::DrawMeshTri3Selective_FaceNorm(
          vtx_xyz, tri_vtx,
          [&tri_flg](unsigned int itri){ return tri_flg[itri] == 0; });
        ::glEnable(GL_LIGHTING);
        ::glColor3d(0,0,0);
        delfem2::opengl::DrawMeshTri3D_Edge(vtx_xyz, tri_vtx);
        ::glColor3d(1, 0, 0);
        Draw(vtx_xyz, tri_vtx, samples);
        glfwSwapBuffers(viewer.window);
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
      if (glfwWindowShouldClose(viewer.window)) { break; }
    }
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
