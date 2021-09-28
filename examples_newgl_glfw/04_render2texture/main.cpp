/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#else
#  include <glad/glad.h>
#endif
#if defined(_MSC_VER)
#  include <windows.h>
#endif
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/shdr_mshtex.h"
#include "delfem2/opengl/new/shdr_mshtri.h"
#include "delfem2/noise.h"
#include "delfem2/mshprimitive.h"

namespace dfm2 = delfem2;

// ---------------------------
dfm2::opengl::CShader_TriMesh shdr0;
dfm2::opengl::CShader_MeshTex shdr_mshtex;
delfem2::glfw::CViewer2 viewer;
unsigned int idTexColor = 0;
unsigned int idTexDepth = 0;

// ---------------------------

void draw(GLFWwindow *window) {
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
//  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  glEnable(GL_TEXTURE_2D);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

  int nw, nh;
  glfwGetFramebufferSize(window, &nw, &nh);
  const float asp = (float) nw / (float)nh;

  {
    glActiveTexture(GL_TEXTURE0); // activate the texture unit first before binding texture
    glBindTexture(GL_TEXTURE_2D, idTexColor);
    delfem2::CMat4f mP, mMV;
    viewer.Mat4_ModelView_Projection(mMV.data(), mP.data(), asp);
    delfem2::CMat4f mZ = delfem2::CMat4f::ScaleXYZ(1.f,1.f,-1.f);
    mMV(0,3) -= 0.5;
    shdr_mshtex.Draw((mZ*mP).transpose().data(), mMV.transpose().data());
  }

  {
    glActiveTexture(GL_TEXTURE0); // activate the texture unit first before binding texture
    glBindTexture(GL_TEXTURE_2D, idTexDepth);
    delfem2::CMat4f mP, mMV;
    viewer.Mat4_ModelView_Projection(mMV.data(), mP.data(), asp);
    delfem2::CMat4f mZ = delfem2::CMat4f::ScaleXYZ(1.f,1.f,-1.f);
    mMV(0,3) += 0.5;
    shdr_mshtex.Draw((mZ*mP).transpose().data(), mMV.transpose().data());
  }

  viewer.SwapBuffers();
  glfwPollEvents();
}

int main() {
  dfm2::glfw::InitGLNew();
  viewer.InitGL();

#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

  {
    std::vector<double> vtx_xyz;
    std::vector<unsigned int> tri_vtx;
    dfm2::MeshTri3_Torus(vtx_xyz, tri_vtx, 0.5, 0.5, 10, 12);
    shdr0.Compile();
    shdr0.Initialize(vtx_xyz, 3, tri_vtx);
  }

  shdr_mshtex.Compile();
  {
    std::vector<float> vtx_xy = {
        -0.5, -0.5,
        +0.5, -0.5,
        +0.5, +0.5,
        -0.5, +0.5,
    };
    std::vector<unsigned int> tri_vtx = {
        0, 1, 2,
        0, 2, 3,
    };
    std::vector<float> vtx_uv = {
        0.0, 0.0,
        1.0, 0.0,
        1.0, 1.0,
        0.0, 1.0
    };
    shdr_mshtex.setCoords(vtx_xy, 2);
    shdr_mshtex.setTexCoords(vtx_uv);
    shdr_mshtex.setElement(tri_vtx, GL_TRIANGLES);
  }

  const unsigned int targetTextureWidth = 256;
  const unsigned int targetTextureHeight = 256;

  {
    ::glEnable(GL_TEXTURE_2D);
    ::glActiveTexture(GL_TEXTURE0);
    // create to render to
    ::glGenTextures(1, &idTexColor);
    ::glBindTexture(GL_TEXTURE_2D, idTexColor);
    // define size and format of level 0
    ::glTexImage2D(GL_TEXTURE_2D, 0,
                   GL_RGBA, targetTextureWidth, targetTextureHeight, 0,
                   GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    // set the filtering so we don't need mips
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }

  {
    const int nSize = 256;
    std::vector<double> aV;
    dfm2::ComputePerlin(aV,
                        nSize, nSize,
                        4, 4, 0.8);
    assert(aV.size() == nSize * nSize);
    std::vector<unsigned char> image(nSize * nSize * 3);
    for (unsigned int i = 0; i < aV.size(); i++) {
      int ival = int(aV[i] * 50 + 128);
      if (ival < 0) { ival = 0; }
      if (ival >= 256) { ival = 255; }
      image[i * 3 + 0] = ival;
      image[i * 3 + 1] = ival;
      image[i * 3 + 2] = ival;
    }
    // -------
    glBindTexture(GL_TEXTURE_2D, idTexColor);
    glTexImage2D(GL_TEXTURE_2D, 0,
                 GL_RGB, targetTextureWidth, targetTextureHeight, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, image.data());
  }

  { // depth texture
    ::glEnable(GL_TEXTURE_2D);
    ::glActiveTexture(GL_TEXTURE0);
    // create to render to
    ::glGenTextures(1, &idTexDepth);
    ::glBindTexture(GL_TEXTURE_2D, idTexDepth);
    // define size and format of level 0
    ::glTexImage2D(GL_TEXTURE_2D, 0,
                   GL_DEPTH_COMPONENT32F, targetTextureWidth, targetTextureHeight, 0,
                   GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    // set the filtering so we don't need mips
//    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }

  {
    // Create and bind the framebuffer
    unsigned int fb;
    ::glGenFramebuffers(1, &fb);
    ::glBindFramebuffer(GL_FRAMEBUFFER, fb);

    // attach the texture as the first color attachment
    ::glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                             idTexColor, 0);
    ::glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
                             idTexDepth, 0);

    // Always check that our framebuffer is ok
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE) {
      std::cout << "error!: " << status << std::endl;
      std::cout << GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT << std::endl;
      std::cout << GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT << std::endl;
      std::cout << GL_FRAMEBUFFER_UNSUPPORTED << std::endl;
      std::cout << GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER << std::endl;
      std::cout << GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE << std::endl;
      std::cout << GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT << std::endl;
      std::cout << GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER << std::endl;
      std::cout << GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER << std::endl;
      return 0;
    }

    int viewport[4];
    ::glGetIntegerv(GL_VIEWPORT, viewport);
    ::glClearColor(0.8, 1.0, 1.0, 1.0);
    ::glClear(GL_DEPTH_BUFFER_BIT);
    ::glViewport(0, 0, targetTextureWidth, targetTextureHeight);
//    ::glClearColor(0.8, 0.8, 0.8, 1.0);
//    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    float mP[16] = {1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1};
    float mMV[16] = {1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1};
    //::glDisable(GL_CULL_FACE);
    ::glEnable(GL_DEPTH_TEST);
    shdr0.Draw(mP, mMV);
    ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
    ::glViewport(0, 0, viewport[2], viewport[3]);
  }

  /*
  {
   std::vector<float> aDepth;
   aDepth.resize(targetTextureHeight*targetTextureWidth);
    ::glBindTexture(GL_TEXTURE_2D, idTexDepth);

    glGetTexImage(GL_TEXTURE_2D,
        0,
        GL_DEPTH_COMPONENT, GL_FLOAT,
                  aDepth.data());
    for(int i=0;i<targetTextureHeight*targetTextureWidth;++i){
      std::cout << i << " " << aDepth[i] << std::endl;
    }
  }
   */

#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

