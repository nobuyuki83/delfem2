/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <glad/glad.h>
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/funcs.h"
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>

// ----------------------------------------------

std::string LoadFile(
    const std::string& fname)
{
  std::ifstream inputFile1(fname.c_str());
  std::istreambuf_iterator<char> vdataBegin(inputFile1);
  std::istreambuf_iterator<char> vdataEnd;
  return std::string(vdataBegin, vdataEnd);
}

void DrawRectangle_FullCanvas()
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glBegin(GL_QUADS);
  glVertex2d(-1, -1);
  glVertex2d(+1, -1);
  glVertex2d(+1, +1);
  glVertex2d(-1, +1);
  glEnd();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glPopAttrib();
}

int main(int argc,char* argv[])
{
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  viewer.camera.view_height = 2.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  std::cout<<"Vendor:"<<glGetString(GL_VENDOR)<<std::endl;
  std::cout<<"GPU: "<<glGetString(GL_RENDERER)<<std::endl;
  std::cout<<"OpenGL ver.: "<<glGetString(GL_VERSION)<<std::endl;
  std::cout<<"OpenGL shading ver.: " <<glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;
  glfwSetWindowTitle(viewer.window, "hoge");

  int id_shader0;
  {
    std::string glslVert = LoadFile(std::string(PATH_INPUT_DIR) + "/glsl.vert");
    std::string glslFrag = LoadFile(std::string(PATH_INPUT_DIR) + "/glsl0.frag");
    id_shader0 = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
  }

  int id_shader1;
  {
    std::string glslVert = LoadFile(std::string(PATH_INPUT_DIR) + "/glsl.vert");
    std::string glslFrag = LoadFile(std::string(PATH_INPUT_DIR) + "/glsl1.frag");
    id_shader1 = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
  }

  while ( !glfwWindowShouldClose(viewer.window) )
  {
    ::glUseProgram(id_shader0);
    for(unsigned int iframe=0;iframe<100;++iframe){
      if( glfwWindowShouldClose(viewer.window) ){ break; }
      GLint iloc = glGetUniformLocation(id_shader0, "resolution");
      GLint viewport[4];
      ::glGetIntegerv(GL_VIEWPORT, viewport);
      glUniform2f(iloc, (float) viewport[2], (float) viewport[3]);
      viewer.DrawBegin_oldGL();
      DrawRectangle_FullCanvas();
      viewer.SwapBuffers();
      glfwPollEvents();
    }
    ::glUseProgram(id_shader1);
    for(unsigned int iframe=0;iframe<100;++iframe){
      if( glfwWindowShouldClose(viewer.window) ){ break; }
      GLint iloc = glGetUniformLocation(id_shader1, "resolution");
      GLint viewport[4];
      ::glGetIntegerv(GL_VIEWPORT, viewport);
      glUniform2f(iloc, (float) viewport[2], (float) viewport[3]);
      iloc = glGetUniformLocation(id_shader1, "time");
      glUniform1f(iloc, glfwGetTime() );
      viewer.DrawBegin_oldGL();
      DrawRectangle_FullCanvas();
      viewer.SwapBuffers();
      glfwPollEvents();
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


