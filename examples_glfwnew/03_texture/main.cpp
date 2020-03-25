/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include "delfem2/imgio.h"
#include "delfem2/noise.h"

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/opengl/gl_funcs.h"
#include "delfem2/opengl/gl_tex.h"
#include "delfem2/opengl/glnew_mshcolor.h"
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glfw_cam.h"

namespace dfm2 = delfem2;

// ---------------------------
CShader_TriMesh_Tex shdr;
delfem2::opengl::CViewer_GLFW viewer;
GLuint m_texName = -1;

// ---------------------------

void draw(GLFWwindow* window)
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//  ::glEnable(GL_DEPTH_TEST);
//  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );

  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0); // activate the texture unit first before binding texture
  glBindTexture(GL_TEXTURE_2D , m_texName);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );

  float mP[16], mMV[16];
  viewer.nav.Matrix_MVP(mMV, mP, window);
  shdr.Draw(mP, mMV);
  
  viewer.DrawEnd_oldGL();
}

int main()
{
  viewer.Init_newGL();
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  shdr.Compile();
  {
    std::vector<double> aPos3d = {
        -0.5, -0.5, 0.0,
        +0.5, -0.5, 0.0,
        +0.5, +0.5, 0.0,
        -0.5, +0.5, 0.0
    };
    std::vector<unsigned int> aTri = {
        0,1,2,
        0,2,3,
    };
    std::vector<double> aTex2d = {
        0.0, 0.0,
        1.0, 0.0,
        1.0, 1.0,
        0.0, 1.0
    };
    shdr.Initialize(aPos3d, aTri, aTex2d);
  }

  {
   const int nSize = 256;
    std::vector<double> aV;
    dfm2::ComputePerlin(aV,
                        nSize, nSize,
                        4, 4, 0.8);
    assert( aV.size() == nSize*nSize );
    std::vector<unsigned char> image(nSize*nSize*3);
    for(unsigned int i=0;i<aV.size();i++){
      int ival = int( aV[i]*50 + 128 );
      if( ival < 0 ){ ival = 0; }
      if( ival >= 256 ){ ival = 255; }
      image[i*3+0] = ival;
      image[i*3+1] = ival;
      image[i*3+2] = ival;
    }
    // -------------
    glEnable(GL_TEXTURE_2D);
    glActiveTexture(GL_TEXTURE0);
    m_texName = dfm2::opengl::SetTexture_RGB(nSize,nSize,image);
    glGenerateMipmap(GL_TEXTURE_2D);
  }

  viewer.nav.camera.view_height = 1.0;
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

