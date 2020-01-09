/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <math.h>
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mshmisc.h"
// ------

#include "glad/glad.h"
#if defined(__APPLE__) && defined(__MACH__)
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include "delfem2/opengl/gl_framebuffer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_color.h"
#include "../glut_cam.h"

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

CNav3D_GLUT nav;
CGLBuffer glbuff;
CElemBuffObj ebo_face;
CElemBuffObj ebo_edge;

// -------------------------------

void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  nav.SetGL_Camera();
  
  delfem2::opengl::DrawBackground(delfem2::CColor(0.8, 0.8, 1.0));
  
  ::glEnable(GL_LIGHTING);
  glbuff.Draw_Start();
  ebo_face.DrawBuffer();
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ebo_edge.DrawBuffer();
  glbuff.Draw_End();
  
  ShowFPS();
  ::glutSwapBuffers();
}

void myGlutIdle(){
  
  if( is_animation ){
  }
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  ::glViewport(0,0,w,h);
  ::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  nav.glutMotion(x, y);
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button, state, x, y);
  ::glutPostRedisplay();
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'a':
      is_animation = !is_animation;
      break;
      
    case '1':
      break;
    case '2':
      break;
    case '3':
      break;
    case '4':
      break;
    case 'i': // one iteration
      break;
    case 'd': // change draw mode
      break;
    case 'f': //
      break;
    case 's': //
      break;
    case ' ':
    {
      static int ifile = 0;
      ifile++;
      if( ifile >= 8 ){ ifile=0; }
    }
  }
  ::glutPostRedisplay();
}


int main(int argc,char* argv[])
{
  glutInit(&argc, argv);
  
  // Initialize GLUT window 3D
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_STENCIL);
  glutCreateWindow("3D View");
  glutDisplayFunc(myGlutDisplay);
  glutIdleFunc(myGlutIdle);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  
  // ---------------------
  
  nav.camera.view_height = 1.0;
  nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Obj(std::string(PATH_INPUT_DIR)+"/bunny_1k.obj",
                    aXYZ, aTri);
  delfem2::Normalize_Points3D(aXYZ);
  std::vector<double> aNorm(aXYZ.size());
  delfem2::Normal_MeshTri3D(aNorm.data(),
                            aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  std::vector<unsigned int> aLine;
  MeshLine_MeshElem(aLine,
                    aTri.data(), aTri.size()/3, delfem2::MESHELEM_TRI, aXYZ.size()/3);
  
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }

  glbuff.SetBuffer_Vtx(aXYZ,3);
  glbuff.SetBuffer_Nrm(aNorm);
  
  ebo_face.SetBuffer_Elem(aTri, GL_TRIANGLES);
  ebo_edge.SetBuffer_Elem(aLine, GL_LINES);
  
  delfem2::opengl::setSomeLighting();
  
  glutMainLoop();
  return 0;
}


