/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <array>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cagedef.h"
#include "delfem2/gizmo_geo3.h"
#include "delfem2/msh_iomisc.h"
#include "delfem2/points.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
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

std::vector<double> ComputeMeanValueCoordinate(
    const std::vector<double>& vec_xyz_ini,
    const std::vector<double>& vec_xyz_cage0,
    const std::vector<unsigned int>& vec_tri_cage){
  const unsigned int num_point = static_cast<unsigned int>(vec_xyz_ini.size() / 3);
  const unsigned int num_point_cage0 = static_cast<unsigned int>(vec_xyz_cage0.size() / 3);
  std::vector<double> matrix0(num_point * num_point_cage0, 0.0);
  for (unsigned int iq = 0; iq < num_point; ++iq) {
    dfm2::CVec3d q0(vec_xyz_ini.data() + iq * 3);
    for (unsigned int itc = 0; itc < vec_tri_cage.size() / 3; ++itc) {
      const unsigned int ip0 = vec_tri_cage[itc * 3 + 0];
      const unsigned int ip1 = vec_tri_cage[itc * 3 + 1];
      const unsigned int ip2 = vec_tri_cage[itc * 3 + 2];
      dfm2::CVec3d p0 = dfm2::CVec3d(vec_xyz_cage0.data() + ip0 * 3) - q0;
      dfm2::CVec3d p1 = dfm2::CVec3d(vec_xyz_cage0.data() + ip1 * 3) - q0;
      dfm2::CVec3d p2 = dfm2::CVec3d(vec_xyz_cage0.data() + ip2 * 3) - q0;
      double w[3];
      dfm2::MeanValueCoordinate_Triangle<dfm2::CVec3d>(w, p0, p1, p2);
      matrix0[iq * num_point_cage0 + ip0] += w[0];
      matrix0[iq * num_point_cage0 + ip1] += w[1];
      matrix0[iq * num_point_cage0 + ip2] += w[2];
    }
  }
  return matrix0;
}

void Example1(
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri) {

  const auto[vec_xyz_cage, vec_tri_cage] =
  MeshTri3_Cuboid(
      {-0.6, -0.6, -0.6},
      {+0.6, +0.6, +0.6});

  const auto num_point = static_cast<unsigned int>(aXYZ0.size() / 3);
  const auto num_point_cage = static_cast<unsigned int>(vec_xyz_cage.size() / 3);
  std::vector<double> matrix = ComputeMeanValueCoordinate(aXYZ0, vec_xyz_cage, vec_tri_cage);
  for (unsigned int iq = 0; iq < num_point; ++iq) {
    double sum_val = 0.0;
    for (unsigned int ip = 0; ip < num_point_cage; ++ip) {
      sum_val += matrix[iq * num_point_cage + ip];
    }
    for (unsigned int ip = 0; ip < num_point_cage; ++ip) {
      matrix[iq * num_point_cage + ip] /= sum_val;
    }
  }
  // --------------------
  delfem2::glfw::CViewer3 viewer;
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::Quat_Bryant(viewer.camera.Quat_tball, -M_PI * 0.25, 0., 0.);
  //
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  std::vector<double> aXYZ = aXYZ0;
  std::vector<double> aXYZ_cage_def = vec_xyz_cage;
  // --------------------
  while (true) {
    const double time = glfwGetTime();
    for (unsigned int ip = 0; ip < num_point_cage; ++ip) {
      aXYZ_cage_def[ip * 3 + 0] = vec_xyz_cage[ip * 3 + 0] + 0.1 * sin(time * ip);
      aXYZ_cage_def[ip * 3 + 1] = vec_xyz_cage[ip * 3 + 1] + 0.1 * sin(time * ip + M_PI * 2 / 3);
      aXYZ_cage_def[ip * 3 + 2] = vec_xyz_cage[ip * 3 + 2] + 0.1 * sin(time * ip + M_PI * 4 / 3);
    }
    for (unsigned int iq = 0; iq < num_point; ++iq){
      aXYZ[iq * 3 + 0] = 0.;
      aXYZ[iq * 3 + 1] = 0.;
      aXYZ[iq * 3 + 2] = 0.;
      for (unsigned int ip = 0; ip < num_point_cage; ++ip) {
        aXYZ[iq * 3 + 0] += matrix[iq * num_point_cage + ip] * aXYZ_cage_def[ip * 3 + 0];
        aXYZ[iq * 3 + 1] += matrix[iq * num_point_cage + ip] * aXYZ_cage_def[ip * 3 + 1];
        aXYZ[iq * 3 + 2] += matrix[iq * num_point_cage + ip] * aXYZ_cage_def[ip * 3 + 2];
      }
    }
    //
    viewer.DrawBegin_oldGL();
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshTri3D_Edge(
        aXYZ_cage_def.data(), aXYZ_cage_def.size() / 3,
        vec_tri_cage.data(), vec_tri_cage.size() / 3);
    delfem2::opengl::DrawMeshTri3D_Edge(
        aXYZ.data(), aXYZ.size() / 3,
        aTri.data(), aTri.size() / 3);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(
        aXYZ.data(),
        aTri.data(), aTri.size() / 3);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
}

void Example2(
    const std::vector<double> &vec_xyz_ini,
    const std::vector<unsigned int> &vec_tri) {
  const auto[vec_xyz_cage0, aTri_cage0] =
  MeshTri3_Cuboid(
      {-0.3, -0.3, 0.2},
      {+0.3, +0.3, 0.4});
  const auto[vec_xyz_cage1, aTri_cage1] =
  MeshTri3_Cuboid(
      {-0.3, -0.3, -0.4},
      {+0.3, +0.3, -0.2});

  const auto num_point = static_cast<unsigned int>(vec_xyz_ini.size() / 3);
  const unsigned int num_point_cage0 = static_cast<unsigned int>(vec_xyz_cage0.size() / 3);
  const unsigned int num_point_cage1 = static_cast<unsigned int>(vec_xyz_cage0.size() / 3);
  std::vector<double> matrix0 = ComputeMeanValueCoordinate(vec_xyz_ini, vec_xyz_cage0, aTri_cage0);
  std::vector<double> matrix1 = ComputeMeanValueCoordinate(vec_xyz_ini, vec_xyz_cage1, aTri_cage1);
  //
  for (unsigned int iq = 0; iq < num_point; ++iq) {
    double sum_val = 0.0;
    for (unsigned int ip = 0; ip < num_point_cage0; ++ip) {
      sum_val += matrix0[iq * num_point_cage0 + ip];
    }
    for (unsigned int ip = 0; ip < num_point_cage1; ++ip) {
      sum_val += matrix1[iq * num_point_cage1 + ip];
    }
    for (unsigned int ip = 0; ip < num_point_cage0; ++ip) {
      matrix0[iq * num_point_cage0 + ip] /= sum_val;
    }
    for (unsigned int ip = 0; ip < num_point_cage1; ++ip) {
      matrix1[iq * num_point_cage1 + ip] /= sum_val;
    }
  }
  // --------------------
  delfem2::glfw::CViewer3 viewer;
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::Quat_Bryant(viewer.camera.Quat_tball, -M_PI * 0.25, 0., 0.);
  //
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  std::vector<double> aXYZ = vec_xyz_ini;
  std::vector<double> vec_xyz_cage0_def = vec_xyz_cage0;
  // --------------------
  while (true) {
    const double time = glfwGetTime();
    for (unsigned int ip = 0; ip < num_point_cage0; ++ip) {
      vec_xyz_cage0_def[ip * 3 + 0] = vec_xyz_cage0[ip * 3 + 0] + 0.1 * sin(time);
      vec_xyz_cage0_def[ip * 3 + 1] = vec_xyz_cage0[ip * 3 + 1] + 0.1 * sin(time * 2);
      vec_xyz_cage0_def[ip * 3 + 2] = vec_xyz_cage0[ip * 3 + 2] + 0.1 * sin(time * 3);
    }
    for (unsigned int iq = 0; iq < num_point; ++iq) {
      aXYZ[iq * 3 + 0] = 0.;
      aXYZ[iq * 3 + 1] = 0.;
      aXYZ[iq * 3 + 2] = 0.;
      for (unsigned int ip = 0; ip < num_point_cage0; ++ip) {
        aXYZ[iq * 3 + 0] += matrix0[iq * num_point_cage0 + ip] * vec_xyz_cage0_def[ip * 3 + 0];
        aXYZ[iq * 3 + 1] += matrix0[iq * num_point_cage0 + ip] * vec_xyz_cage0_def[ip * 3 + 1];
        aXYZ[iq * 3 + 2] += matrix0[iq * num_point_cage0 + ip] * vec_xyz_cage0_def[ip * 3 + 2];
      }
      for (unsigned int ip = 0; ip < num_point_cage1; ++ip) {
        aXYZ[iq * 3 + 0] += matrix1[iq * num_point_cage1 + ip] * vec_xyz_cage1[ip * 3 + 0];
        aXYZ[iq * 3 + 1] += matrix1[iq * num_point_cage1 + ip] * vec_xyz_cage1[ip * 3 + 1];
        aXYZ[iq * 3 + 2] += matrix1[iq * num_point_cage1 + ip] * vec_xyz_cage1[ip * 3 + 2];
      }
    }
    //
    viewer.DrawBegin_oldGL();
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vec_xyz_cage0_def.data(), vec_xyz_cage0_def.size() / 3,
        aTri_cage0.data(), aTri_cage0.size() / 3);
    delfem2::opengl::DrawMeshTri3D_Edge(
        vec_xyz_cage1.data(), vec_xyz_cage1.size() / 3,
        aTri_cage1.data(), aTri_cage1.size() / 3);
    delfem2::opengl::DrawMeshTri3D_Edge(
        aXYZ.data(), aXYZ.size() / 3,
        vec_tri.data(), vec_tri.size() / 3);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(
        aXYZ.data(),
        vec_tri.data(), vec_tri.size() / 3);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
}


// ---------------

int main() {
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Ply(
      std::string(PATH_INPUT_DIR) + "/bunny_1k.ply",
      aXYZ, aTri);
  delfem2::Normalize_Points3(aXYZ);
  Example1(aXYZ, aTri);
  Example2(aXYZ, aTri);
}


