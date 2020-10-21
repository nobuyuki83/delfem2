/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#include "delfem2/primitive.h"
#include "delfem2/points.h"

// gl related includes
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

#ifndef M_PI
#  define M_PI 3.141592
#endif

namespace dfm2 = delfem2;

// ---------------------------------------------------

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aElm;
  // ----------------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  viewer.nav.camera.Rot_Camera(-0.5,-0.5);
  dfm2::opengl::setSomeLighting();

  while (true)
  {
    // -----
    dfm2::MeshTri3_Capsule(aXYZ,aElm,
        0.2, 1.0,
        16, 5, 8);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aElm);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // --------------
    dfm2::MeshQuad2D_Grid(aXYZ,aElm,16,8);
    dfm2::Scale_Points(aXYZ.data(),aXYZ.size()/3,3,0.1);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshQuad2D_Edge(aXYZ, aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // --------------
    {
      const double pmin[3] = {-0.8, -0.2, -0.5};
      const double pmax[3] = {+0.8, +0.2, +0.5};
      dfm2::MeshQuad3_CubeVox(aXYZ, aElm, pmin,pmax );
    }
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshQuad3D_FaceNorm(aXYZ.data(), aElm.data(), aElm.size()/4);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshQuad3D_Edge(aXYZ, aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // --------------
    dfm2::MeshTri3_Torus(aXYZ,aElm,
        0.5,0.1,32,32);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aElm);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // ----
    dfm2::MeshTri3D_Cube(aXYZ,aElm,8);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aElm);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // ----
    dfm2::MeshTri3D_Sphere(aXYZ,aElm,
        0.5, 32,32);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aElm);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // ----
    dfm2::MeshTri3D_Icosahedron(aXYZ,aElm);
    dfm2::Scale_Points(aXYZ.data(),
        aXYZ.size()/3,3,0.5);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aElm);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // ----
    dfm2::MeshTri3D_CylinderOpen(aXYZ,aElm,
        0.2,1.0, 16, 8);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aElm);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // -----
    dfm2::MeshTri3D_CylinderClosed(aXYZ,aElm,
        0.2,1.0, 16, 8);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aElm);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
    // -----
    dfm2::MeshTri3D_Disk(aXYZ,aElm,
        0.5, 8, 32);
    for(int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aElm);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aElm);
      viewer.DrawEnd_oldGL();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
