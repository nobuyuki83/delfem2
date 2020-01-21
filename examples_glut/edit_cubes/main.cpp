/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <complex>
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/voxel.h"
#include "delfem2/camera.h"

// -------

#if defined(__APPLE__) && defined(__MACH__)
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_color.h"
#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glold_voxbv.h"
#include "../glut_cam.h"

#ifndef M_PI
  #define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// -----------------------------------------

CNav3D_GLUT nav;
int imode_draw = 0;

enum EDIT_MODE {
  EDIT_NONE,
  EDIT_ADD,
  EDIT_DELETE,
};
EDIT_MODE edit_mode = EDIT_NONE;

int imode_sym = 2;

const double elen = 1.0;
const CVector3 org(0,0,0);
std::vector<CCubeGrid> aCubeGrid;
int icube_picked;
int iface_picked;

// -----------------------------------------

void myGlutDisplay(void)
{
  ::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  //	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  delfem2::opengl::DrawBackground(delfem2::CColor(0.2, 0.7, 0.7));
  
  nav.SetGL_Camera();
  
  CVector3 offsym(0,0,0);
  if( imode_sym == 2 ){ offsym.p[2] = -elen*0.5; }
  for(int ic=0;ic<aCubeGrid.size();++ic){
    Draw_CubeGrid(ic==icube_picked, iface_picked, elen, org+offsym, aCubeGrid[ic]);
  }
  
  if( edit_mode == EDIT_ADD ){
    DrawMessage_LeftTop("mode add");
  }
  if( edit_mode == EDIT_DELETE ){
    DrawMessage_LeftTop("mode delete");
  }
  
  
  ::glLineWidth(1);
  delfem2::opengl::DrawAxis(1);
  
  ShowFPS();
  
  ::glutSwapBuffers();
}

void myGlutIdle(){
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  ::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key,x,y);
}

void myGlutMotion( int x, int y ){
  nav.glutMotion(x,y);
  if( nav.imodifier != 0){ return; }
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  dfm2::CVector2 sp0(nav.mouse_x-nav.dx, nav.mouse_y-nav.dy);
  dfm2::CVector2 sp1(nav.mouse_x, nav.mouse_y);
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button,state,x,y);
  if( nav.imodifier != 0 ) return;
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  dfm2::CVector2 sp0(nav.mouse_x, nav.mouse_y);
  const CVector3 src_pick = screenUnProjection(CVector3(sp0.x(),sp0.y(), 0.0), mMV,mPj);
  const CVector3 dir_pick = screenUnProjectionDirection(CVector3(0.0,  0, -1.0 ), mMV,mPj);
  if( state == GLUT_DOWN ){
    CVector3 offsym(0,0,0);
    if( imode_sym == 2 ){ offsym.p[2] = -elen*0.5; }
    double src_pick0[3]; src_pick.CopyValueTo(src_pick0);
    double dir_pick0[3]; dir_pick.CopyValueTo(dir_pick0);
    double offsym0[3];   offsym.CopyValueTo(offsym0);
    Pick_CubeGrid(icube_picked, iface_picked,
                  src_pick0,dir_pick0, elen, offsym0, aCubeGrid);
    if( edit_mode == EDIT_ADD ){
      int ivx1,ivy1,ivz1;
      Adj_CubeGrid(ivx1,ivy1,ivz1,
                   icube_picked,iface_picked,aCubeGrid);
      Add_CubeGrid(aCubeGrid,ivx1,ivy1,ivz1);
      if( imode_sym == 1 ){ Add_CubeGrid(aCubeGrid,ivx1,ivy1,-ivz1-1); }
      if( imode_sym == 2 ){ Add_CubeGrid(aCubeGrid,ivx1,ivy1,-ivz1); }
    }
    if( edit_mode == EDIT_DELETE ){
      int ivx1 = aCubeGrid[icube_picked].ivx;
      int ivy1 = aCubeGrid[icube_picked].ivy;
      int ivz1 = aCubeGrid[icube_picked].ivz;
      Del_CubeGrid(aCubeGrid,ivx1,ivy1,ivz1);
      if( imode_sym == 1 ){ Del_CubeGrid(aCubeGrid,ivx1,ivy1,-ivz1-1); }
      if( imode_sym == 2 ){ Del_CubeGrid(aCubeGrid,ivx1,ivy1,-ivz1); }
    }
  }
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
      break;
    case '0':
    {
      edit_mode = EDIT_NONE;
      std::cout << "edit_none" << std::endl;
      break;
    }
    case '1':
    {
      edit_mode = EDIT_ADD;
      std::cout << "edit_add" << std::endl;
      break;
    }
    case '2':
    {
      edit_mode = EDIT_DELETE;
      std::cout << "edit_delete" << std::endl;
      break;
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
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("3D View");
  glutDisplayFunc(myGlutDisplay);
  glutIdleFunc(myGlutIdle);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  
  // -------------------------
  nav.camera.view_height = 2.0;
  nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
   
  delfem2::opengl::setSomeLighting();
  
  aCubeGrid.push_back( CCubeGrid(0,0,0) );
  
  glutMainLoop();
  return 0;
}


