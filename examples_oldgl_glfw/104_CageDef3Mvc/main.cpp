/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cagedef.h"
#include "delfem2/gizmo_geo3.h"
#include "delfem2/mshio.h"
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

int main(
    [[maybe_unused]] int argc,
    [[maybe_unused]] char *argv[]) {
  delfem2::glfw::CViewer3 viewer;
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Ply(
      std::string(PATH_INPUT_DIR) + "/bunny_1k.ply",
      aXYZ, aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<unsigned int> aTri_cage;
  std::vector<double> aXYZ_cage;
  {
    double pmin[3], pmax[3];
    delfem2::BoundingBox3_Points3(
        pmin, pmax,
        aXYZ.data(), aXYZ.size() / 3);
    pmin[0] -= 0.1;
    pmin[1] -= 0.1;
    pmin[2] -= 0.1;
    pmax[0] += 0.1;
    pmax[1] += 0.1;
    pmax[2] += 0.1;
    std::vector<unsigned int> aQuad_cage;
    delfem2::MeshQuad3_CubeVox(
        aXYZ_cage, aQuad_cage,
        pmin, pmax);
    delfem2::convert2Tri_Quad(
        aTri_cage, aQuad_cage);
  }
  const auto num_point = static_cast<unsigned int>(aXYZ.size() / 3);
  const auto num_point_cage = static_cast<unsigned int>(aXYZ_cage.size() / 3);
  std::vector<double> matrix(num_point * num_point_cage, 0.0);
  for (unsigned int iq = 0; iq < num_point; ++iq) {
    dfm2::CVec3d q0(aXYZ.data() + iq * 3);
    for (unsigned int itc = 0; itc < aTri_cage.size() / 3; ++itc) {
      unsigned int ip0 =  aTri_cage[itc * 3 + 0];
      unsigned int ip1 =  aTri_cage[itc * 3 + 1];
      unsigned int ip2 =  aTri_cage[itc * 3 + 2];
      dfm2::CVec3d p0 = dfm2::CVec3d(aXYZ_cage.data() + ip0 * 3) - q0;
      dfm2::CVec3d p1 = dfm2::CVec3d(aXYZ_cage.data() + ip1 * 3) - q0;
      dfm2::CVec3d p2 = dfm2::CVec3d(aXYZ_cage.data() + ip2 * 3) - q0;
      double w[3];
      dfm2::MeanValueCoordinate_Triangle<dfm2::CVec3d>(w, p0, p1, p2);
      matrix[iq * num_point_cage + ip0] += w[0];
      matrix[iq * num_point_cage + ip1] += w[1];
      matrix[iq * num_point_cage + ip2] += w[2];
    }
    double sum_val = 0.0;
    for(unsigned int ip=0;ip<num_point_cage;++ip) {
      sum_val += matrix[iq * num_point_cage + ip];
    }
    for(unsigned int ip=0;ip<num_point_cage;++ip) {
      matrix[iq * num_point_cage + ip] /= sum_val;
    }
  }
    // --------------------
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::Quat_Bryant(viewer.camera.Quat_tball, -M_PI*0.25, 0., 0.);
  delfem2::opengl::setSomeLighting();
  const std::vector<double> aXYZ0_cage = aXYZ_cage;
  // --------------------
  while (true) {
    const double time = glfwGetTime();
    for(unsigned int ip=0;ip<num_point_cage;++ip){
      aXYZ_cage[ip*3+0] = aXYZ0_cage[ip*3+0] + 0.1*sin(time*ip);
      aXYZ_cage[ip*3+1] = aXYZ0_cage[ip*3+1] + 0.1*sin(time*ip+M_PI*2/3);
      aXYZ_cage[ip*3+2] = aXYZ0_cage[ip*3+2] + 0.1*sin(time*ip+M_PI*4/3);
    }
    for (unsigned int iq = 0; iq < num_point; ++iq) {
      double p[3] = {0,0,0};
      for(unsigned int ip=0;ip<num_point_cage;++ip){
        p[0] += matrix[iq*num_point_cage + ip] * aXYZ_cage[ip*3+0];
        p[1] += matrix[iq*num_point_cage + ip] * aXYZ_cage[ip*3+1];
        p[2] += matrix[iq*num_point_cage + ip] * aXYZ_cage[ip*3+2];
      }
      aXYZ[iq*3+0] = p[0];
      aXYZ[iq*3+1] = p[1];
      aXYZ[iq*3+2] = p[2];
    }
    //
    viewer.DrawBegin_oldGL();
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshTri3D_Edge(
        aXYZ_cage.data(), aXYZ_cage.size() / 3,
        aTri_cage.data(), aTri_cage.size() / 3);
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
  exit(EXIT_SUCCESS);
}


