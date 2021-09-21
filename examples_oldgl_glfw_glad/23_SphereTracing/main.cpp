/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/funcs.h"

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

int main()
{
  delfem2::glfw::CViewer3 viewer;
  viewer.projection
  = std::make_unique<delfem2::Projection_LookOriginFromZplus<double>>(1.0, true, 45.0);
  //
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  delfem2::opengl::setSomeLighting();
  std::cout<<"Vendor:"<<glGetString(GL_VENDOR)<<std::endl;
  std::cout<<"GPU: "<<glGetString(GL_RENDERER)<<std::endl;
  std::cout<<"OpenGL ver.: "<<glGetString(GL_VERSION)<<std::endl;
  std::cout<<"OpenGL shading ver.: " <<glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;
  glfwSetWindowTitle(viewer.window, "hoge");

  int id_shader0;
  {
    std::string glslVert = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl.vert");
    std::string glslFrag = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl0.frag");
    id_shader0 = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
  }

  int id_shader1;
  {
    std::string glslVert = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl.vert");
    std::string glslFrag = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl1.frag");
    id_shader1 = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
  }

  int id_shader2;
  {
    std::string glslVert = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl.vert");
    std::string glslFrag = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl2.frag");
    id_shader2 = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
  }

  while ( !glfwWindowShouldClose(viewer.window) )
  {
    ::glUseProgram(id_shader0);
    for(unsigned int iframe=0;iframe<100;++iframe){
      if( glfwWindowShouldClose(viewer.window) ){ break; }
      GLint iloc = ::glGetUniformLocation(id_shader0, "resolution");
      GLint viewport[4];
      ::glGetIntegerv(GL_VIEWPORT, viewport);
      ::glUniform2f(iloc, (float) viewport[2], (float) viewport[3]);
      viewer.DrawBegin_oldGL();
      DrawRectangle_FullCanvas();
      viewer.SwapBuffers();
      glfwPollEvents();
    }
    // ---------------
    ::glUseProgram(id_shader1);
    for(unsigned int iframe=0;iframe<100;++iframe){
      if( glfwWindowShouldClose(viewer.window) ){ break; }
      GLint iloc = glGetUniformLocation(id_shader1, "resolution");
      GLint viewport[4];
      ::glGetIntegerv(GL_VIEWPORT, viewport);
      glUniform2f(iloc, (float) viewport[2], (float) viewport[3]);

      iloc = glGetUniformLocation(id_shader2, "mMVPinv");
      float mMV[16], mMVP[16], mMVPinv[16];
      delfem2::CMat4f mP = viewer.projection->Mat4ColumnMajor((float) viewport[2] / (float) viewport[3]);
      viewer.modelview.Mat4ColumnMajor(mMV);
      delfem2::MatMat4(mMVP,mMV,mP.data());
      delfem2::Inverse_Mat4(mMVPinv,mMVP);
      glUniformMatrix4fv(iloc,1,GL_FALSE,mMVPinv);
      iloc = glGetUniformLocation(id_shader2, "mMV");
      glUniformMatrix4fv(iloc,1,GL_FALSE,mMV);
      viewer.DrawBegin_oldGL();
      DrawRectangle_FullCanvas();
      viewer.SwapBuffers();
      glfwPollEvents();
    }
    // ---------
    ::glUseProgram(id_shader2);
    for(unsigned int iframe=0;iframe<100;++iframe){
      if( glfwWindowShouldClose(viewer.window) ){ break; }
      GLint iloc = glGetUniformLocation(id_shader2, "resolution");
      GLint viewport[4];
      ::glGetIntegerv(GL_VIEWPORT, viewport);
      glUniform2f(iloc, (float) viewport[2], (float) viewport[3]);
      iloc = glGetUniformLocation(id_shader2, "mMVPinv");
      float mMV[16], mMVP[16], mMVPinv[16];
      delfem2::CMat4f mP = viewer.projection->Mat4ColumnMajor((float) viewport[2] / (float) viewport[3]);
      viewer.modelview.Mat4ColumnMajor(mMV);
      delfem2::MatMat4(mMVP,mMV,mP.data());
      delfem2::Inverse_Mat4(mMVPinv,mMVP);
      glUniformMatrix4fv(iloc,1,GL_FALSE,mMVPinv);
      iloc = glGetUniformLocation(id_shader2, "mMV");
      glUniformMatrix4fv(iloc,1,GL_FALSE,mMV);
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


