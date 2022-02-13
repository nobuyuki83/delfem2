/*
 * Copyright (c) 2020-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/msh_primitive.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/openglstb/glyph.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

int main()
{
  delfem2::openglstb::CGlyph glyph(
      std::string(PATH_INPUT_DIR)+"/myFont.png");
  glyph.ParseGlyphInfo(
      std::string(PATH_INPUT_DIR)+"/myFont.fnt");
  // -------------------
  
  dfm2::glfw::CViewer3 viewer(2.0);
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  std::cout<<"Vendor:"<<glGetString(GL_VENDOR)<<'\n';
  std::cout<<"GPU: "<<glGetString(GL_RENDERER)<<'\n';
  std::cout<<"OpenGL ver.: "<<glGetString(GL_VERSION)<<'\n';
  // GL_SHADING_LANGUAGE_VERSION is undefined on Windows
  // std::cout<<"OpenGL shading ver.: " <<glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;

  glyph.InitGL();
  
  // ----------------------
  while (!glfwWindowShouldClose(viewer.window))
  {
    glfwSetWindowTitle(viewer.window, "naive");
    //
    for(int iframe=0;iframe<100;++iframe){
      viewer.DrawBegin_oldGL();
      glyph.DrawEntireGlyph();
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
    }
    //
    for(char ichar=30;ichar<90;++ichar) {
      for (int iframe=0; iframe<5; ++iframe) {
        viewer.DrawBegin_oldGL();
        glyph.DrawCharAt(ichar, 0.01, 0, 0);
        glfwSwapBuffers(viewer.window);
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) goto EXIT;
      }
    }
    //
    for (int iframe=0; iframe<100; ++iframe) {
      viewer.DrawBegin_oldGL();
      glyph.DrawStringAt("Hello DelFEM2!", 0.01, -1, 0);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if (glfwWindowShouldClose(viewer.window)) goto EXIT;
    }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


