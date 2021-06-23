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
#include "delfem2/opengl/tex.h"
#include "delfem2/imgio.h"
#include "delfem2/mshprimitive.h"
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <cmath>


namespace dfm2 = delfem2;

// ----------------------------

// -----------------------------

std::string LoadFile(
    const std::string& fname)
{
  std::ifstream inputFile1(fname.c_str());
  std::istreambuf_iterator<char> vdataBegin(inputFile1);
  std::istreambuf_iterator<char> vdataEnd;
  return std::string(vdataBegin,vdataEnd);
}

void setShaderProgram(
    int& id_shader_program,
    std::string glslVert,
    std::string glslFrag)
{
  glUseProgram(0);
  glDeleteProgram(id_shader_program);
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
  
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  viewer.camera.view_height = 2.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
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
  
  
  while (true)
  {
    int id_shapder_program;
    /*
    {
      glfwSetWindowTitle(viewer.window, "phong");
      const std::string glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_adsphong.vert");
      const std::string glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_adsphong.frag");
      setShaderProgram(id_shapder_program, glslVert,glslFrag);
      for(int iframe=0;iframe<100;++iframe){
        viewer.DrawBegin_oldGL();
        myGlutDisplay(id_shapder_program);
        viewer.DrawEnd_oldGL();
        if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
      }
    }
    {
      glfwSetWindowTitle(viewer.window, "gouraud");
      const std::string glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_adsgouraud.vert");
      const std::string glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_adsgouraud.frag");
      setShaderProgram(id_shapder_program, glslVert,glslFrag);
      for(int iframe=0;iframe<100;++iframe){
        viewer.DrawBegin_oldGL();
        myGlutDisplay(id_shapder_program);
        viewer.DrawEnd_oldGL();
        if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
      }
    }
    {
      glfwSetWindowTitle(viewer.window, "toon");
      const std::string glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_toon.vert");
      const std::string glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_toon.frag");
      setShaderProgram(id_shapder_program, glslVert,glslFrag);
      for(int iframe=0;iframe<100;++iframe){
        viewer.DrawBegin_oldGL();
        myGlutDisplay(id_shapder_program);
        viewer.DrawEnd_oldGL();
        if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
      }
    }
     */
    {
      glfwSetWindowTitle(viewer.window, "normalmap");
      const std::string glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_normalmap.vert");
      const std::string glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_normalmap.frag");
      setShaderProgram(id_shapder_program, glslVert,glslFrag);
      for(int iframe=0;iframe<100;++iframe){
        viewer.DrawBegin_oldGL();
        myGlutDisplay(id_shapder_program);
        viewer.SwapBuffers();
        glfwPollEvents();
        if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
      }
    }
    /*
    {
      glfwSetWindowTitle(viewer.window, "texture");
      const std::string glslVert = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_simpletexture.vert");
      const std::string glslFrag = LoadFile(std::string(PATH_INPUT_DIR)+"/glsl120_simpletexture.frag");
      setShaderProgram(id_shapder_program, glslVert,glslFrag);
      for(int iframe=0;iframe<100;++iframe){
        viewer.DrawBegin_oldGL();
        myGlutDisplay(id_shapder_program);
        viewer.DrawEnd_oldGL();
        if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
      }
    }
     */
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


