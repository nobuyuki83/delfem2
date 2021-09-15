/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "delfem2/msh_io_obj.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

// -------------------------------

class CElemBuffObj {
 public:
  void SetBuffer_Elem(
      const std::vector<unsigned int> &aTri,
      unsigned int gl_elem_type);
  void DrawBuffer() const;
 public:
  unsigned int iebo;
  unsigned int gl_elem_type;
  unsigned int size_elem;
  bool is_lighting;
};

void CElemBuffObj::SetBuffer_Elem(
    const std::vector<unsigned int> &aTri,
    unsigned int gl_elem_type_) {
  this->gl_elem_type = gl_elem_type_;
  size_elem = aTri.size();
  glGenBuffers(1, &iebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               aTri.size() * (sizeof(int)),
               aTri.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void CElemBuffObj::DrawBuffer() const {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->iebo);
  glDrawElements(gl_elem_type, size_elem, GL_UNSIGNED_INT, nullptr);
}

// ----------

class CGLBuffer {
 public:
  CGLBuffer() {
    vbo = -1;
    vbo_nrm = -1;
    ndim = 0;
  }
  void SetBuffer_Vtx(const std::vector<double> &aXYZ, int ndim);
  void SetBuffer_Nrm(const std::vector<double> &aNrm);
  void Draw_Start() const;
  static void Draw_End();
 public:

 public:
  unsigned int vbo;
  unsigned int vbo_nrm;
  unsigned int ndim;
};

void CGLBuffer::SetBuffer_Vtx(
    const std::vector<double> &aXYZ,
    int ndim_) {
  this->ndim = ndim_;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER,
               aXYZ.size() * (sizeof(double)),
               aXYZ.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void CGLBuffer::SetBuffer_Nrm(const std::vector<double> &aNrm) {
  glGenBuffers(1, &vbo_nrm);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_nrm);
  glBufferData(GL_ARRAY_BUFFER,
               aNrm.size() * (sizeof(double)),
               aNrm.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void CGLBuffer::Draw_Start() const {
  assert(glIsBuffer(this->vbo));
  glEnableClientState(GL_VERTEX_ARRAY); // this makes glold
  glBindBuffer(GL_ARRAY_BUFFER, this->vbo);
  glVertexPointer(ndim, GL_DOUBLE, 0, 0);
  if (glIsBuffer(this->vbo_nrm)) {
    glEnableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_nrm);
    glNormalPointer(GL_DOUBLE, 0, 0);
  }
}
void CGLBuffer::Draw_End() {
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);
}



// -------------------------------

bool is_animation;

CGLBuffer glbuff;
CElemBuffObj ebo_face;
CElemBuffObj ebo_edge;

// -------------------------------

void myGlutDisplay() {
  delfem2::opengl::DrawBackground(delfem2::CColor(0.8, 0.8, 1.0));
  ::glEnable(GL_LIGHTING);
  glbuff.Draw_Start();
  ebo_face.DrawBuffer();
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0, 0, 0);
  ebo_edge.DrawBuffer();
  glbuff.Draw_End();
}

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  delfem2::Read_Obj(
      vtx_xyz, tri_vtx,
      std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.obj");
  delfem2::Normalize_Points3(vtx_xyz);
  std::vector<double> vtx_norm(vtx_xyz.size());
  delfem2::Normal_MeshTri3D(
      vtx_norm.data(),
      vtx_xyz.data(), vtx_xyz.size() / 3,
      tri_vtx.data(), tri_vtx.size() / 3);
  std::vector<unsigned int> aLine;
  MeshLine_MeshElem(
      aLine,
      tri_vtx.data(), tri_vtx.size() / 3, delfem2::MESHELEM_TRI,
      vtx_xyz.size() / 3);

  // ---------------------------

  dfm2::glfw::CViewer3 viewer;
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  viewer.camera.Rot_Camera(-0.5, -0.5);
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  if (!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  dfm2::opengl::setSomeLighting();

  glbuff.SetBuffer_Vtx(vtx_xyz, 3);
  glbuff.SetBuffer_Nrm(vtx_norm);

  ebo_face.SetBuffer_Elem(tri_vtx, GL_TRIANGLES);
  ebo_edge.SetBuffer_Elem(aLine, GL_LINES);

  while (true) {
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
  }

  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


