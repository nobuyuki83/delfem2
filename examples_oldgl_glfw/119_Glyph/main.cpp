/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/openglstb/glyph.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/str.h"
#include <iostream>
#include <fstream>
#include <cmath>

namespace dfm2 = delfem2;

// ------------------------------------------------------


// ------------------------------------------------------

int main(int argc,char* argv[])
{
  delfem2::openglstb::CGlyph glyph(std::string(PATH_INPUT_DIR)+"/myFont.png");
  glyph.ParseGlyphInfo(std::string(PATH_INPUT_DIR)+"/myFont.fnt");
  // -------------------
  
  dfm2::opengl::CViewer3 viewer;
  viewer.Init_oldGL();
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  viewer.camera.view_height = 2.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  std::cout<<"Vendor:"<<glGetString(GL_VENDOR)<<'\n';
  std::cout<<"GPU: "<<glGetString(GL_RENDERER)<<'\n';
  std::cout<<"OpenGL ver.: "<<glGetString(GL_VERSION)<<'\n';
  std::cout<<"OpenGL shading ver.: " <<glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;

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


