/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mshmisc.h"
// ------
#include "glad/glad.h"
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"
namespace dfm2 = delfem2;

// -------------------------------

class CElemBuffObj{
public:
  void SetBuffer_Elem(const std::vector<unsigned int>& aTri, unsigned int gl_elem_type);
  void DrawBuffer() const ;
public:
  unsigned int iebo;
  unsigned int gl_elem_type;
  unsigned int size_elem;
  bool is_lighting;
};


void CElemBuffObj::SetBuffer_Elem(const std::vector<unsigned int>& aTri, unsigned int gl_elem_type)
{
  this->gl_elem_type = gl_elem_type;
  size_elem = aTri.size();
  glGenBuffers(1,&iebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               aTri.size()*(sizeof(int)),
               aTri.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void CElemBuffObj::DrawBuffer() const
{
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->iebo);
  glDrawElements(gl_elem_type, size_elem, GL_UNSIGNED_INT, 0);
}

// ----------

class CGLBuffer
{
public:
  CGLBuffer(){
    vbo = -1;
    vbo_nrm = -1;
    ndim = 0;
  }
  void SetBuffer_Vtx(const std::vector<double>& aXYZ, int ndim);
  void SetBuffer_Nrm(const std::vector<double>& aNrm);
  void Draw_Start() const;
  void Draw_End() const ;
public:
  
public:
  unsigned int vbo;
  unsigned int vbo_nrm;
  unsigned int ndim;
};


void CGLBuffer::SetBuffer_Vtx(const std::vector<double>& aXYZ, int ndim)
{
  this->ndim = ndim;
  glGenBuffers(1,&vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER,
               aXYZ.size() * (sizeof(double)),
               aXYZ.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void CGLBuffer::SetBuffer_Nrm(const std::vector<double>& aNrm)
{
  glGenBuffers(1,&vbo_nrm);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_nrm);
  glBufferData(GL_ARRAY_BUFFER,
               aNrm.size() * (sizeof(double)),
               aNrm.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void CGLBuffer::Draw_Start() const
{
  assert( glIsBuffer(this->vbo) );
  glEnableClientState(GL_VERTEX_ARRAY); // this makes glold
  glBindBuffer(GL_ARRAY_BUFFER, this->vbo);
  glVertexPointer(ndim, GL_DOUBLE, 0, 0);
  if( glIsBuffer(this->vbo_nrm) ){
    glEnableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_nrm);
    glNormalPointer(GL_DOUBLE, 0, 0);
  }
}
void CGLBuffer::Draw_End() const
{
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);
}



// -------------------------------

bool is_animation;

CGLBuffer glbuff;
CElemBuffObj ebo_face;
CElemBuffObj ebo_edge;

// -------------------------------

void myGlutDisplay()
{
  delfem2::opengl::DrawBackground(delfem2::CColor(0.8, 0.8, 1.0));
  ::glEnable(GL_LIGHTING);
  glbuff.Draw_Start();
  ebo_face.DrawBuffer();
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ebo_edge.DrawBuffer();
  glbuff.Draw_End();
}

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Obj(std::string(PATH_INPUT_DIR)+"/bunny_1k.obj",
                    aXYZ, aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<double> aNorm(aXYZ.size());
  delfem2::Normal_MeshTri3D(aNorm.data(),
                            aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  std::vector<unsigned int> aLine;
  MeshLine_MeshElem(aLine,
                    aTri.data(), aTri.size()/3, delfem2::MESHELEM_TRI, aXYZ.size()/3);

  // ---------------------------

  dfm2::opengl::CViewer_GLFW viewer;
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  viewer.nav.camera.Rot_Camera(-0.5,-0.5);
  viewer.Init_oldGL();
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  dfm2::opengl::setSomeLighting();

  glbuff.SetBuffer_Vtx(aXYZ,3);
  glbuff.SetBuffer_Nrm(aNorm);
  
  ebo_face.SetBuffer_Elem(aTri, GL_TRIANGLES);
  ebo_edge.SetBuffer_Elem(aLine, GL_LINES);

  while(true){
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }

EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


