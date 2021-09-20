/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <glad/glad.h>
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/shdr_msh.h"
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
  return {vdataBegin, vdataEnd};
}

int main()
{
  delfem2::glfw::CViewer3 viewer;
  viewer.projection.view_height = 1.0;
  viewer.projection.is_pars = true;
  viewer.projection.fovy = 45.0;
  delfem2::opengl::CShader_Mesh shdr;
  shdr.color[0] = 1;
  //
  delfem2::glfw::InitGLNew();
  viewer.InitGL();
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
  //  delfem2::opengl::setSomeLighting();
  std::cout<<"Vendor:"<<glGetString(GL_VENDOR)<<std::endl;
  std::cout<<"GPU: "<<glGetString(GL_RENDERER)<<std::endl;
  std::cout<<"OpenGL ver.: "<<glGetString(GL_VERSION)<<std::endl;
  std::cout<<"OpenGL shading ver.: " <<glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;
  glfwSetWindowTitle(viewer.window, "Sphere Tracing Demo");

//  shdr.Compile();
  {
    std::vector<float> aXY = {
        -1, -1,
        +1, -1,
        +1, +1,
        -1, +1,
    };
    std::vector<unsigned int> aTri0 = {
        0,1,2,
        0,2,3,
    };
    shdr.Initialize(aXY, 2, aTri0, GL_TRIANGLES );
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
    int id_shader = -1;
    if( (int)(glfwGetTime()*0.5) % 2 == 0 ) {
      id_shader = id_shader2;
    } else{
      id_shader = id_shader1;
    }
    //
    ::glUseProgram(id_shader);
    {
      if( glfwWindowShouldClose(viewer.window) ){ break; }
      GLint iloc = glGetUniformLocation(id_shader, "resolution");
      GLint viewport[4];
      ::glGetIntegerv(GL_VIEWPORT, viewport);
      glUniform2f(iloc, (float) viewport[2], (float) viewport[3]);
      iloc = glGetUniformLocation(id_shader, "mMVPinv");
      float mMV[16], mP[16], mMVP[16], mMVPinv[16];
      viewer.projection.Mat4ColumnMajor(mP, (float) viewport[2] / (float) viewport[3]);
      viewer.modelview.Mat4ColumnMajor(mMV);
      delfem2::MatMat4(mMVP,mMV,mP);
      delfem2::Inverse_Mat4(mMVPinv,mMVP);
      glUniformMatrix4fv(iloc,1,GL_FALSE,mMVPinv);
      iloc = glGetUniformLocation(id_shader, "mMV");
      glUniformMatrix4fv(iloc,1,GL_FALSE,mMV);
      shdr.vao.Draw(0);
      viewer.SwapBuffers();
      glfwPollEvents();
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


