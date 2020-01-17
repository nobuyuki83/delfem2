/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include "delfem2/noise.h"
#include "delfem2/primitive.h"
#include "delfem2/vec3.h"

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/opengl/glnew_mshcolor.h"
#include "delfem2/opengl/glnew_gpusampler.h"
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/glfw_cam.h"


namespace dfm2 = delfem2;

// ---------------------------
CShader_TriMesh shdr0;
CShader_Points shdr1;
delfem2::opengl::CViewer_GLFW viewer;
CGPUSampler sampler;

// ---------------------------

void draw(GLFWwindow* window)

{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );

  float mP[16], mMV[16];
  viewer.nav.Matrix_MVP(mMV, mP, window);
  shdr1.Draw(mP,mMV);
  shdr0.Draw(mP, mMV);
  sampler.Draw(mP,mMV);
  
  viewer.DrawEnd_oldGL();
}

int main()
{
  {
    int nres = 256;
    double elen = 0.01;
    sampler.Init(nres, nres, "4byte",false);
    sampler.SetCoord(elen, 4.0,
                     CVector3(-nres*elen*0.5,nres*elen*0.5,-2).stlvec(),
                     CVector3(0,0,-1).stlvec(),
                     CVector3(1,0,0).stlvec() );
    sampler.SetColor(1, 0, 0);
    sampler.draw_len_axis = 1.0;
    sampler.bgcolor = {1,1,1};
  }

  viewer.Init_newGL();
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  sampler.InitGL();

  {
    std::vector<double> aXYZ;
    std::vector<unsigned int> aTri;
    dfm2::MeshTri3_Torus(aXYZ, aTri, 0.8, 0.1, 8, 8);
    shdr0.Compile();
    shdr0.Initialize(aXYZ, aTri);
  }

  sampler.Start();
  float mMV[16], mP[16];
  sampler.Matrix_MVP(mMV,mP);
  shdr0.Draw(mP,mMV);
  sampler.End();

  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

