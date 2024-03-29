/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should be before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "delfem2/vec3.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_normal.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/opengl/funcs.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/r2tglo.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

namespace dfm2 = delfem2;

// -----------------------------

std::string LoadFile(
    const std::string &fname) {
  std::ifstream inputFile1(fname.c_str());
  std::istreambuf_iterator<char> vdataBegin(inputFile1);
  std::istreambuf_iterator<char> vdataEnd;
  return {vdataBegin, vdataEnd};
}

// ------------------------------------------------------

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  dfm2::Read_Obj(
      vtx_xyz, tri_vtx,
      std::filesystem::path(PATH_ASSET_DIR) / "bunny_1k.obj");
  dfm2::Normalize_Points3(vtx_xyz, 2.5);
  dfm2::Rotate_Points3(
      vtx_xyz,
      -M_PI * 0.5, 0.0, 0.0);
  std::vector<double> vtx_normal(vtx_xyz.size());
  dfm2::Normal_MeshTri3D(
      vtx_normal.data(),
      vtx_xyz.data(), vtx_xyz.size() / 3,
      tri_vtx.data(), tri_vtx.size() / 3);
  // ---------------------------------------

  unsigned int nres = 256;
  double elen = 0.02;
  dfm2::opengl::CRender2Tex sampler;
  sampler.SetTextureProperty(nres, nres, true);
  dfm2::CMat4d::AffineAxisTransform(
      {1,0,0},
      {0,1,0},
      {0,0,1}
  ).CopyTo(sampler.mat_modelview);
  dfm2::CMat4d::AffineOrthogonalProjection(
      -elen*nres*0.5, elen*nres*0.5,
      -elen*nres*0.5, elen*nres*0.5,
      -elen*nres*0.5, elen*nres*0.5
  ).CopyTo(sampler.mat_projection);
  dfm2::opengl::CDrawerOldGL_Render2Tex draw_sampler;
  draw_sampler.SetPointColor(1, 0, 0);
  draw_sampler.draw_len_axis = 1.0;

  // -------------------
  dfm2::glfw::CViewer3 viewer(4.0);
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  if (!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  delfem2::opengl::setSomeLighting();
  std::cout << "Vendor:" << glGetString(GL_VENDOR) << std::endl;
  std::cout << "GPU: " << glGetString(GL_RENDERER) << std::endl;
  std::cout << "OpenGL ver.: " << glGetString(GL_VERSION) << std::endl;
  std::cout << "OpenGL shading ver.: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
  glfwSetWindowTitle(viewer.window, "naive");

  // -------------------------

  unsigned int id_tex_depth;
  ::glGenTextures(1, &id_tex_depth);
  ::glBindTexture(GL_TEXTURE_2D, id_tex_depth);

  int id_shader_normal;
  {
    std::string glslVert = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl120_normal.vert");
    std::string glslFrag = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl120_normal.frag");
    id_shader_normal = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
    glUseProgram(0);
  }

  /*
  int id_shader_laplace;
  {
    std::string glslVert = LoadFile(std::string(PATH_INPUT_DIR) + "/glsl120_normal.vert");
    std::string glslFrag = LoadFile(std::string(PATH_INPUT_DIR) + "/glsl120_normal.frag");
    id_shader_normal = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
    glUseProgram(0);
  }
   */
  ::glActiveTexture(0);
  sampler.InitGL();
  unsigned int id_tex_norm;
  { // normal information to norm tex
    sampler.Start();
    dfm2::opengl::SetView(sampler);
    ::glClearColor(0.0, 0.0, 0.0, 1.0);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glUseProgram(id_shader_normal);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz, tri_vtx, vtx_normal);
    sampler.End();
    // ---------
    ::glGenTextures(1, &id_tex_norm);
    ::glActiveTexture(1);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_norm);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    // ---------
    ::glTexImage2D(GL_TEXTURE_2D,
                   0, GL_RGBA,
                   static_cast<int>(sampler.width),
                   static_cast<int>(sampler.height),
                   0, GL_RGBA,
                   GL_UNSIGNED_BYTE, sampler.aRGBA_8ui.data());
  }

  { // rendering image normally
    sampler.Start();
    ::glClearColor(1.0, 1.0, 1.0, 1.0);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glUseProgram(0);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz, tri_vtx);
    sampler.End();
  }

  int id_shader_edge;
  {
    std::string glslVert = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl120_edge.vert");
    std::string glslFrag = LoadFile(std::string(PATH_SOURCE_DIR) + "/glsl120_edge.frag");
    id_shader_edge = delfem2::opengl::setUpGLSL(glslVert, glslFrag);
    ::glUseProgram(id_shader_edge);
    {
      GLint texLoc0 = glGetUniformLocation(id_shader_edge, "TexOrg");
      glUniform1i(texLoc0, 0); // GL_TEXTURE0
    }
    {
      GLint texLoc0 = glGetUniformLocation(id_shader_edge, "TexNrm");
      glUniform1i(texLoc0, 1); // GL_TEXTURE0
    }
    {
      GLint texLoc0 = glGetUniformLocation(id_shader_edge, "nTexWidth");
      glUniform1i(texLoc0, static_cast<int>(sampler.width)); // GL_TEXTURE0
    }
    {
      GLint texLoc0 = glGetUniformLocation(id_shader_edge, "nTexHeight");
      glUniform1i(texLoc0, static_cast<int>(sampler.height)); // GL_TEXTURE0
    }
    glUseProgram(0);
  }

  // -----
  draw_sampler.isDrawTex = true;

  while (true) {
    viewer.DrawBegin_oldGL();
    { // draw mesh
      ::glEnable(GL_LIGHTING);
      glUseProgram(0);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz, tri_vtx);
    }
    { // draw texture
      ::glUseProgram(id_shader_edge);
      ::glActiveTexture(GL_TEXTURE0);
      ::glBindTexture(GL_TEXTURE_2D, sampler.id_tex_color);
      ::glActiveTexture(GL_TEXTURE1);
      ::glBindTexture(GL_TEXTURE_2D, id_tex_norm);
      draw_sampler.Draw_Texture(sampler);
    }
    { // draw bounding box
      ::glUseProgram(0);
      ::glDisable(GL_LIGHTING);
      ::glLineWidth(1);
      ::glColor3d(0, 0, 0);
      draw_sampler.Draw_BoundingBox(sampler);
      draw_sampler.Draw_Axis(sampler);
    }
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) goto EXIT;
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


