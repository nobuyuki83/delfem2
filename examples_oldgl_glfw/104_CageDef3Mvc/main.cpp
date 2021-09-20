/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <array>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cagedef.h"
#include "delfem2/gizmo_geo3.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/points.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/mshuni.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/vec3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

#ifndef M_PI
#  define M_PI 3.141592
#endif

namespace dfm2 = delfem2;

// -------------------------------------------

std::pair<std::vector<double>, std::vector<unsigned int> >
MeshTri3_Cuboid(
    const std::array<double, 3> &pmin,
    const std::array<double, 3> &pmax) {
  std::vector<double> vec_xyz;
  std::vector<unsigned int> vec_quad_cage;
  delfem2::MeshQuad3_CubeVox(
      vec_xyz, vec_quad_cage,
      pmin.data(), pmax.data());
  std::vector<unsigned int> vec_tri;
  delfem2::convert2Tri_Quad(
      vec_tri, vec_quad_cage);
  return {vec_xyz, vec_tri};
}

void Example1(
    const std::vector<double> &vtx_xyz0,
    const std::vector<unsigned int> &tri_vtxidx) {

  const auto[vec_xyz_cage, vec_tri_cage] =
  MeshTri3_Cuboid(
      {-0.6, -0.6, -0.6},
      {+0.6, +0.6, +0.6});

  const size_t num_vtx = vtx_xyz0.size() / 3;
  const size_t num_vtx_cage = vec_xyz_cage.size() / 3;
  std::vector<double> matrix = delfem2::ComputeMeanValueCoordinate<dfm2::CVec3d>(
      vtx_xyz0, vec_xyz_cage, vec_tri_cage);
  for (unsigned int iq = 0; iq < num_vtx; ++iq) {
    double sum_val = 0.0;
    for (unsigned int ip = 0; ip < num_vtx_cage; ++ip) {
      sum_val += matrix[iq * num_vtx_cage + ip];
    }
    for (unsigned int ip = 0; ip < num_vtx_cage; ++ip) {
      matrix[iq * num_vtx_cage + ip] /= sum_val;
    }
  }
  // --------------------
  delfem2::glfw::CViewer3 viewer;
  viewer.projection.view_height = 1.0;
  delfem2::Quat_Bryant(viewer.modelview.Quat_tball, -M_PI * 0.25, 0., 0.);
  //
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  std::vector<double> vtx_xyz = vtx_xyz0;
  std::vector<double> vtx_xyz_cage_def = vec_xyz_cage;
  // --------------------
  while (true) {
    const double time = glfwGetTime();
    for (unsigned int ip = 0; ip < num_vtx_cage; ++ip) {
      vtx_xyz_cage_def[ip * 3 + 0] = vec_xyz_cage[ip * 3 + 0] + 0.1 * sin(time * ip);
      vtx_xyz_cage_def[ip * 3 + 1] = vec_xyz_cage[ip * 3 + 1] + 0.1 * sin(time * ip + M_PI * 2 / 3);
      vtx_xyz_cage_def[ip * 3 + 2] = vec_xyz_cage[ip * 3 + 2] + 0.1 * sin(time * ip + M_PI * 4 / 3);
    }
    for (unsigned int iq = 0; iq < num_vtx; ++iq) {
      vtx_xyz[iq * 3 + 0] = 0.;
      vtx_xyz[iq * 3 + 1] = 0.;
      vtx_xyz[iq * 3 + 2] = 0.;
      for (unsigned int ip = 0; ip < num_vtx_cage; ++ip) {
        vtx_xyz[iq * 3 + 0] += matrix[iq * num_vtx_cage + ip] * vtx_xyz_cage_def[ip * 3 + 0];
        vtx_xyz[iq * 3 + 1] += matrix[iq * num_vtx_cage + ip] * vtx_xyz_cage_def[ip * 3 + 1];
        vtx_xyz[iq * 3 + 2] += matrix[iq * num_vtx_cage + ip] * vtx_xyz_cage_def[ip * 3 + 2];
      }
    }
    //
    viewer.DrawBegin_oldGL();
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vtx_xyz_cage_def.data(), vtx_xyz_cage_def.size() / 3,
        vec_tri_cage.data(), vec_tri_cage.size() / 3);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vtx_xyz.data(), vtx_xyz.size() / 3,
        tri_vtxidx.data(), tri_vtxidx.size() / 3);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(
        vtx_xyz.data(),
        tri_vtxidx.data(), tri_vtxidx.size() / 3);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
}

void Example2(
    const std::vector<double> &vtx_xyz_ini,
    const std::vector<unsigned int> &tri_vtx_ind) {
  const auto[vtx_xyz_cage0, tri_vtx_ind_cage0] =
  MeshTri3_Cuboid(
      {-0.3, -0.3, 0.2},
      {+0.3, +0.3, 0.4});
  const auto[vtx_xyz_cage1, tri_vtx_ind_cage1] =
  MeshTri3_Cuboid(
      {-0.3, -0.3, -0.4},
      {+0.3, +0.3, -0.2});

  std::vector<double> weight_dof_dofcage{1, 1, 1, 1, 1, 1, 1, 1};

  const size_t num_vtx = vtx_xyz_ini.size() / 3;
  std::vector<double> vector0 = delfem2::ComputeMeanValueCoordinateReduced<dfm2::CVec3d>(
      vtx_xyz_ini, vtx_xyz_cage0, tri_vtx_ind_cage0,
      1, weight_dof_dofcage);
  std::vector<double> vector1 = delfem2::ComputeMeanValueCoordinateReduced<dfm2::CVec3d>(
      vtx_xyz_ini, vtx_xyz_cage1, tri_vtx_ind_cage1,
      1, weight_dof_dofcage);
  for (unsigned int iq = 0; iq < num_vtx; ++iq) {
    double sum_val = vector0[iq] + vector1[iq];
    vector0[iq] /= sum_val;
    vector1[iq] /= sum_val;
  }
  for (unsigned int iq = 0; iq < num_vtx; ++iq) {
    double sum_val = vector0[iq] + vector1[iq];
    vector0[iq] /= sum_val;
    vector1[iq] /= sum_val;
  }
  // --------------------
  delfem2::glfw::CViewer3 viewer;
  viewer.projection.view_height = 1.0;
  delfem2::Quat_Bryant(viewer.modelview.Quat_tball, -M_PI * 0.25, 0., 0.);
  //
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  std::vector<double> vec_xyz;
  std::vector<double> vec_xyz_cage0_def = vtx_xyz_cage0;
  // --------------------
  while (true) {
    const double time = glfwGetTime();
    double disp[3] = {
        0.1 * sin(time),
        0.1 * sin(time * 2),
        0.1 * sin(time * 3)};
    const size_t num_vtx_cage0 = vtx_xyz_cage0.size() / 3;
    for (unsigned int ip = 0; ip < num_vtx_cage0; ++ip) {
      vec_xyz_cage0_def[ip * 3 + 0] = vtx_xyz_cage0[ip * 3 + 0] + disp[0];
      vec_xyz_cage0_def[ip * 3 + 1] = vtx_xyz_cage0[ip * 3 + 1] + disp[1];
      vec_xyz_cage0_def[ip * 3 + 2] = vtx_xyz_cage0[ip * 3 + 2] + disp[2];
    }
    vec_xyz = vtx_xyz_ini;
    for (unsigned int iq = 0; iq < num_vtx; ++iq) {
      vec_xyz[iq * 3 + 0] += vector0[iq] * disp[0];
      vec_xyz[iq * 3 + 1] += vector0[iq] * disp[1];
      vec_xyz[iq * 3 + 2] += vector0[iq] * disp[2];
    }
    //
    viewer.DrawBegin_oldGL();
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vec_xyz_cage0_def.data(), vec_xyz_cage0_def.size() / 3,
        tri_vtx_ind_cage0.data(), tri_vtx_ind_cage0.size() / 3);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vtx_xyz_cage1.data(), vtx_xyz_cage1.size() / 3,
        tri_vtx_ind_cage1.data(), tri_vtx_ind_cage1.size() / 3);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vec_xyz.data(), vec_xyz.size() / 3,
        tri_vtx_ind.data(), tri_vtx_ind.size() / 3);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(
        vec_xyz.data(),
        tri_vtx_ind.data(), tri_vtx_ind.size() / 3);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
}


// ---------------

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  delfem2::Read_Ply(
      vtx_xyz, tri_vtx,
      std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.ply");
  delfem2::Normalize_Points3(vtx_xyz);
  Example1(vtx_xyz, tri_vtx);
  Example2(vtx_xyz, tri_vtx);
}


