/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "delfem2/imgio.h"
#include "delfem2/primitive.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/gl_funcs.h"
#include "delfem2/opengl/tex_gl.h"
#include "delfem2/opengl/funcs_glold.h"

namespace dfm2 = delfem2;

// ----------------------------

// -----------------------------

std::string LoadFile
(const std::string& fname)
{
  std::ifstream inputFile1(fname.c_str());
  std::istreambuf_iterator<char> vdataBegin(inputFile1);
  std::istreambuf_iterator<char> vdataEnd;
  return std::string(vdataBegin,vdataEnd);
}

void setShaderProgram(
    int& id_shader_program,
    int isp)
{
  std::string glslVert, glslFrag;
  glUseProgram(0);
  glDeleteProgram(id_shader_program);
  if( isp == 0 ){
    glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl_adsphong.vert");
    glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl_adsphong.frag");
  }
  else if( isp == 1 ){
    glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl_adsgouraud.vert");
    glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl_adsgouraud.frag");
  }
  else if( isp == 2 ){
    glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl_toon.vert");
    glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl_toon.frag");
  }
  else if( isp == 3 ){
    glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl_simpletexture.vert");
    glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl_simpletexture.frag");
  }
  id_shader_program = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
  //  
  glUseProgram(id_shader_program);
  {
    GLint texLoc = glGetUniformLocation(id_shader_program, "myTexColor");
    glUniform1i(texLoc, 0); // GL_TEXTURE0
  }
  
  {
    GLint texLoc = glGetUniformLocation(id_shader_program, "myTexNormal");
    glUniform1i(texLoc, 1); // GL_TEXTURE1
  }
}

// ------------------------------------------------------

void myGlutDisplay(
    int id_shader_program)
{
  ::glEnable(GL_LIGHTING);
  {
    float lightPosition[4] = { 0.0, 0.0, 5.0, 1.0 };
    float light_ambient[4] = { 0.3, 0.3, 0.3, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  }
  {
    float kd[4] = {1.0, 0.0, 0.0, 1.0};
    float shininess = 100.0;
    float ks[4] = {1.0, 1.0, 1.0, 1.0};
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,kd);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,kd);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS,shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, ks);
  }
  
  glUseProgram(id_shader_program);
  dfm2::opengl::DrawTorus_Solid(1.0, 0.4, 2.0);
  glUseProgram(0);
}


int main(int argc,char* argv[])
{
  
  // -------------------
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  // --------------
  
  std::cout<<"Vendor:"<<glGetString(GL_VENDOR)<<'\n';
  std::cout<<"GPU: "<<glGetString(GL_RENDERER)<<'\n';
  std::cout<<"OpenGL ver.: "<<glGetString(GL_VERSION)<<'\n';
  std::cout<<"OpenGL shading ver.: " <<glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;
  
  // -------------------------
  
  dfm2::SFile_TGA tga_color;  LoadTGAFile(std::string(PATH_INPUT_DIR)+"/rock_color.tga",  &tga_color);
  dfm2::SFile_TGA tga_normal; LoadTGAFile(std::string(PATH_INPUT_DIR)+"/rock_normal.tga", &tga_normal);
  
  GLuint aIndTex[2];
  ::glGenTextures(2, aIndTex);
  
  // color
  ::glActiveTexture(GL_TEXTURE0);
  ::glBindTexture(GL_TEXTURE_2D, aIndTex[0]);
  ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  ::glTexImage2D(GL_TEXTURE_2D,
                 0, GL_RGBA,
                 tga_color.imageWidth, tga_color.imageHeight,
                 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, tga_color.imageData);
  
  // normal
  ::glActiveTexture(GL_TEXTURE1);
  ::glBindTexture(GL_TEXTURE_2D, aIndTex[1]);
  ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  ::glTexImage2D(GL_TEXTURE_2D,
                 0, GL_RGBA,
                 tga_normal.imageWidth, tga_normal.imageHeight,
                 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, tga_normal.imageData);
  
  ::glEnable(GL_TEXTURE_2D);
  
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    int id_shapder_program;
    glfwSetWindowTitle(viewer.window, "phong");
    setShaderProgram(id_shapder_program, 0);
    for(int iframe=0;iframe<100;++iframe){
      viewer.DrawBegin_oldGL();
      myGlutDisplay(id_shapder_program);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
    }
    // ----
    glfwSetWindowTitle(viewer.window, "gouraud");
    setShaderProgram(id_shapder_program, 1);
    for(int iframe=0;iframe<100;++iframe){
      viewer.DrawBegin_oldGL();
      myGlutDisplay(id_shapder_program);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
    }
    // -----
    glfwSetWindowTitle(viewer.window, "toon");
    setShaderProgram(id_shapder_program, 2);
    for(int iframe=0;iframe<100;++iframe){
      viewer.DrawBegin_oldGL();
      myGlutDisplay(id_shapder_program);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
    }
    // -----
    glfwSetWindowTitle(viewer.window, "texture");
    setShaderProgram(id_shapder_program, 3);
    for(int iframe=0;iframe<100;++iframe){
      viewer.DrawBegin_oldGL();
      myGlutDisplay(id_shapder_program);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
    }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


