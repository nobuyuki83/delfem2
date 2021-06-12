/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengleigen/funcs.h"
#include "delfem2/geoconvhull3.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <cstdlib>
#include <random>

namespace dfm2 = delfem2;

// ---------------------------------------

int main(int argc,char* argv[])
{
  std::vector<Eigen::Vector3f,Eigen::aligned_allocator<Eigen::Vector3f> > aXYZ(100);
  std::vector<unsigned int> aTri;
  // --
  delfem2::glfw::CViewer3 viewer;
  viewer.camera.view_height = 1.5;
  // ---
  delfem2::glfw::InitGLOld();
  viewer.InitGL();

  double time_last_update = -2;
  while (!glfwWindowShouldClose(viewer.window))
  {
    const double time_now = glfwGetTime();
    if(time_now - time_last_update > 1 ){
      std::mt19937 rngeng(std::random_device{}());
      std::uniform_real_distribution<double> dist_m1p1(-1,+1);
      for(auto& xyz : aXYZ){
        xyz(0) = dist_m1p1(rngeng);
        xyz(1) = dist_m1p1(rngeng);
        xyz(2) = dist_m1p1(rngeng);
      }
      delfem2::ConvexHull<double>(aTri,aXYZ);
      time_last_update = time_now;
    }
    //
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    ::glEnable(GL_CULL_FACE);
    ::glCullFace(GL_BACK);
    {
      ::glColor3d(0,0,0);
      ::glPointSize(3);
      delfem2::opengleigen::DrawPoints(aXYZ);
      //
      ::glColor4d(1,1,1,0.8);
      delfem2::opengleigen::DrawMeshTri3(aTri,aXYZ);
      //
      ::glColor3d(0,0,0);
      delfem2::opengleigen::DrawMeshTri3_Edge(aTri,aXYZ);
    }
    
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
