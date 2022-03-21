/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <chrono>
#include <vector>
#include <cstdlib>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should put before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/geo_convhull3.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/eigen_opengl/funcs.h"

namespace dfm2 = delfem2;
// ---------------------------------------

template <class VEC, typename REAL = typename VEC::Scalar>
double MeasureTime()
{
  std::vector<VEC,Eigen::aligned_allocator<VEC> > vtx_xyz(1000);
  {
    std::mt19937 rngeng(std::random_device{}());
    std::uniform_real_distribution<REAL> dist_m1p1(-1,+1);
    for(auto& xyz : vtx_xyz){
      xyz.setZero();
      xyz[0] = dist_m1p1(rngeng);
      xyz[1] = dist_m1p1(rngeng);
      xyz[2] = dist_m1p1(rngeng);
    }
  }
  std::vector<unsigned int> tri_vtx;
  auto start = std::chrono::system_clock::now();
  for(int it=0;it<1000;++it) {
    delfem2::ConvexHull3(tri_vtx, vtx_xyz);
  }
  auto end = std::chrono::system_clock::now();
  return std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
}

int main()
{
  std::cout << "Available :SIMD Instructions: "<< Eigen::SimdInstructionSetsInUse() << std::endl;
  std::vector<Eigen::Vector4f,Eigen::aligned_allocator<Eigen::Vector4f> > vtx_xyz(100);
  std::vector<unsigned int> tri_vtx;
  delfem2::glfw::CViewer3 viewer(1.5);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  double time_last_update = -2;
  while (!glfwWindowShouldClose(viewer.window)) {
    const double time_now = glfwGetTime();
    if(time_now - time_last_update > 1 ){
      std::mt19937 rngeng(std::random_device{}());
      std::uniform_real_distribution<float> dist_m1p1(-1,+1);
      for(auto& xyz : vtx_xyz){
        xyz[0] = dist_m1p1(rngeng);
        xyz[1] = dist_m1p1(rngeng);
        xyz[2] = dist_m1p1(rngeng);
      }
      delfem2::ConvexHull3(tri_vtx, vtx_xyz);
      time_last_update = time_now;
    }
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    ::glEnable(GL_CULL_FACE);
    ::glCullFace(GL_BACK);
    {
      ::glColor3d(0,0,0);
      ::glPointSize(3);
      delfem2::eigen_opengl::DrawPoints(vtx_xyz);
      //
      ::glColor4d(1,1,1,0.8);
      delfem2::eigen_opengl::DrawMeshTri3(tri_vtx, vtx_xyz);
      //
      ::glColor3d(0,0,0);
      delfem2::eigen_opengl::DrawMeshTri3_Edge(tri_vtx, vtx_xyz);
    }
    
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  // ----------
  std::cout << MeasureTime<Eigen::Vector3f>() << std::endl;
  std::cout << MeasureTime<Eigen::Vector4f>() << std::endl;
  std::cout << MeasureTime<Eigen::Vector3d>() << std::endl;
  std::cout << MeasureTime<Eigen::Vector4d>() << std::endl;
  exit(EXIT_SUCCESS);
}
