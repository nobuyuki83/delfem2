/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include <math.h>
#include <time.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/msh.h"

#include "delfem2/paramgeo_v23.h"

#include "delfem2/color_gl.h"
#include "delfem2/funcs_glut.h"


CVector3 GetPointSurf
(double u, double v,
 int isurf,
 std::vector<int>& aIndCP,
 std::vector<CVector3>& aCP)
{
  int i00 = aIndCP[isurf*16+ 0];
  int i01 = aIndCP[isurf*16+ 1];
  int i02 = aIndCP[isurf*16+ 2];
  int i03 = aIndCP[isurf*16+ 3];
  int i04 = aIndCP[isurf*16+ 4];
  int i05 = aIndCP[isurf*16+ 5];
  int i06 = aIndCP[isurf*16+ 6];
  int i07 = aIndCP[isurf*16+ 7];
  int i08 = aIndCP[isurf*16+ 8];
  int i09 = aIndCP[isurf*16+ 9];
  int i10 = aIndCP[isurf*16+10];
  int i11 = aIndCP[isurf*16+11];
  int i12 = aIndCP[isurf*16+12];
  int i13 = aIndCP[isurf*16+13];
  int i14 = aIndCP[isurf*16+14];
  int i15 = aIndCP[isurf*16+15];
  return getPointSurfaceBezierCubic(u,v,
                                    aCP[i00], aCP[i01], aCP[i02], aCP[i03],
                                    aCP[i04], aCP[i05], aCP[i06], aCP[i07],
                                    aCP[i08], aCP[i09], aCP[i10], aCP[i11],
                                    aCP[i12], aCP[i13], aCP[i14], aCP[i15]);
}

void AddQuads
(std::vector<CVector3>& aPQuad,
 int n,
 int isurf,
 std::vector<int>& aIndCP,
 std::vector<CVector3>& aCP)
{
  for (int i = 0; i<(n+1); i++){
    for (int j = 0; j<(n+1); j++){
      for(int ic=0;ic<4;++ic){
        int iu = i;
        int jv = j;
        if(ic==1){ iu++; }
        if(ic==2){ iu++; jv++; }
        if(ic==3){ jv++; }
        double u = (double)iu/n;
        double v = (double)jv/n;
        CVector3 p = GetPointSurf(u, v, isurf, aIndCP, aCP);
        aPQuad.push_back(p);
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////

std::vector<int> aIndCP;
std::vector<CVector3> aCP;
std::vector<CVector3> aPQuad;
int n = 20;

bool is_animation = true;
CGlutWindowManager win;

//////////////////////////////////////////////////////////////////////////////////////


static void myGlVertex3d(const CVector3& v)
{
  ::glVertex3d(v.x,v.y,v.z);
}

static void myGlVertex3d(int i, const std::vector<CVector3>& aV)
{
  const CVector3& v = aV[i];
  ::glVertex3d(v.x,v.y,v.z);
}


void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  win.SetGL_Camera();
  
  DrawBackground();
  
  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  ::glEnable(GL_CULL_FACE);
  ::glCullFace(GL_BACK);
  {
    ::glColor3d(0,0,0);
    ::glPointSize(3);
    ::glBegin(GL_POINTS);
    for(int icp=0;icp<(int)aCP.size();icp++){
      myGlVertex3d(aCP[icp]);
    }
    ::glEnd();
  }
  
  {
    int nq = (int)aPQuad.size()/4;
    ::glColor4d(1,1,1,1.0);
    ::glBegin(GL_QUADS);
    for(int iq=0;iq<nq;iq++){
      myGlVertex3d(iq*4+0,aPQuad);
      myGlVertex3d(iq*4+1,aPQuad);
      myGlVertex3d(iq*4+2,aPQuad);
      myGlVertex3d(iq*4+3,aPQuad);
    }
    ::glEnd();
    ////
    ::glColor4d(0,0,0,1.0);
    for(int iq=0;iq<nq;iq++){
      ::glBegin(GL_LINE_STRIP);
      myGlVertex3d(iq*4+0,aPQuad);
      myGlVertex3d(iq*4+1,aPQuad);
      myGlVertex3d(iq*4+2,aPQuad);
      myGlVertex3d(iq*4+3,aPQuad);
      myGlVertex3d(iq*4+0,aPQuad);
      ::glEnd();
    }
  }
  ::glutSwapBuffers();
}


void myGlutIdle()
{
  if( !is_animation ){
    ::glutPostRedisplay();
    return;
  }
  ::glutPostRedisplay();
}

void myGlutResize(int w, int h)
{
  ::glViewport(0, 0, w, h);
  ::glutPostRedisplay();  
}

void myGlutMotion( int x, int y )
{
  win.glutMotion(x,y);
  ::glutPostRedisplay();
}

void myGlutMouse(int ibutton, int state, int x, int y)
{
  win.glutMouse(ibutton, state, x, y);
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
      break;
    case '\033':
      break;
    case 'a':
      is_animation = !is_animation;
      break;
    case ' ':
    {
      int nCP = 16;
      aCP.resize(nCP);
      for(int iCP=0;iCP<16;iCP++){
        aCP[iCP].x = 1.0*(double)rand()/(RAND_MAX+1.0);
        aCP[iCP].y = 1.0*(double)rand()/(RAND_MAX+1.0);
        aCP[iCP].z = 1.0*(double)rand()/(RAND_MAX+1.0);
      }
      aIndCP.resize(16);
      for(int i=0;i<16;++i){ aIndCP[i] = i; }
      aCP[ 0] += CVector3(-2,-2,0);
      aCP[ 1] += CVector3(-2,-1,0);
      aCP[ 2] += CVector3(-2,+1,0);
      aCP[ 3] += CVector3(-2,+2,0);
      ////
      aCP[ 4] += CVector3(-1,-2,0);
      aCP[ 5] += CVector3(-1,-1,0);
      aCP[ 6] += CVector3(-1,+1,0);
      aCP[ 7] += CVector3(-1,+2,0);
      ////
      aCP[ 8] += CVector3(+1,-2,0);
      aCP[ 9] += CVector3(+1,-1,0);
      aCP[10] += CVector3(+1,+1,0);
      aCP[11] += CVector3(+1,+2,0);
      ////
      aCP[12] += CVector3(+2,-2,0);
      aCP[13] += CVector3(+2,-1,0);
      aCP[14] += CVector3(+2,+1,0);
      aCP[15] += CVector3(+2,+2,0);
      
      /////
      aPQuad.clear();
      AddQuads(aPQuad,n,0,aIndCP,aCP);
      break;
    }
  }
  ::glutPostRedisplay();
}


void myGlutSpecial(int Key, int x, int y)
{
  win.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  ::glutInit(&argc, argv);
  
  // Initialize GLUT
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("Initial");
  
  // Define callback functions
  glutDisplayFunc(myGlutDisplay);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  glutIdleFunc(myGlutIdle);
  
  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  win.camera.view_height = 4;
  
  
  glutMainLoop();
  return 0;
}
