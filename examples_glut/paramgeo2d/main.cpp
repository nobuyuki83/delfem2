/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/vec2.h"
#include "delfem2/camera.h"

#include "delfem2/paramgeo_v23.h"

#include "delfem2/gl2_v23.h"
#include "delfem2/gl2_funcs.h"

// ---------------------------------------------

// should be inside cam2
int ipoint_picked;
double pos_down[2];
double pos_cur[2];

int ndegree = 2;
std::vector<int> aKnotMulti;
std::vector<double> aKnot;
std::vector<double> aKnotFlat;
std::vector<CVector2> aCtrlPoint;

const int nsmpl = 100;
std::vector<CVector2> polyline0; // current test

// -----------------------------------------------

void SetExample
(int ndeg, int ncp)
{
  ndegree = ndeg;
//  const int nk = ncp+ndeg+1;
  const int ndiv = ncp-ndeg;
  aKnot.assign(ndiv+1,0);
  for(int idiv=0;idiv<ndiv+1;++idiv){
    aKnot[idiv] = (double)idiv/ndiv;
  }
  aKnotMulti.assign(ndiv+1,1);
  aKnotMulti[   0] = ndeg+1;
  aKnotMulti[ndiv] = ndeg+1;
  FlatKnot(aKnotFlat, aKnotMulti, aKnot);
  for(unsigned int ik=0;ik<aKnotFlat.size();++ik){
    std::cout << "knot" << ik << " " << aKnotFlat[ik] << std::endl;
  }
  /////
  aCtrlPoint.assign(ncp, CVector2(0,0));  //7+2+1 = 10
  for(unsigned int i=0;i<aCtrlPoint.size();++i){
    aCtrlPoint[i].x = i*2.0/(aCtrlPoint.size()-1)-1.0;
  }
}

void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset(1.1f, 4.0f);
  
  setGL_Camera2D();
  
  ::glDisable(GL_LIGHTING);
  
  ::glLineWidth(2);
  ::glPointSize(5);
  ::glColor3d(1,0,0);
  
  ::glColor3d(0,0.5,0);
  drawPolyLine2D(aCtrlPoint);
  
  
  
  ::glPointSize(2);
  ::glColor3d(0,0,0);
  drawPolyLine2D(polyline0);
  
  ::glutSwapBuffers();
}

void myGlutIdle(){
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  ::glViewport(0,0,w,h);
  ::glutPostRedisplay();
}

void myGlutMotion(int x, int y)
{
  getPosOnScreen_Camera2D(pos_cur[0], pos_cur[1],
                          x, y);
  if( ipoint_picked != -1 ){
    aCtrlPoint[ipoint_picked].x = pos_cur[0];
    aCtrlPoint[ipoint_picked].y = pos_cur[1];
  }
  SampleBSpline(polyline0, nsmpl, ndegree, aKnotFlat, aCtrlPoint);
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  if (state==GLUT_DOWN){
    getPosOnScreen_Camera2D(pos_down[0], pos_down[1],
                            x, y);
    ipoint_picked = -1;
    for(unsigned int ip=0;ip<aCtrlPoint.size();ip++){
      if( Distance(CVector2(pos_down[0],pos_down[1]),aCtrlPoint[ip])<0.03 ){
        ipoint_picked = ip;
      }
    }
  }
  else{
    ipoint_picked = -1;
  }
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
      break;
    case 'a': //
      break;
      /*
       case '1':
       {
       --nsmpl_bezier;
       if( nsmpl_bezier == 0 ){
       nsmpl_bezier = 1;
       }
       updatedCV();
       }
       break;
       case '2':
       {
       ++nsmpl_bezier;
       updatedCV();
       }
       break;
       */
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
  glutInitWindowPosition(100,100);
  glutInitWindowSize(1000, 600);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_STENCIL);
  glutCreateWindow("3D View");
  glutDisplayFunc(myGlutDisplay);
  glutIdleFunc(myGlutIdle);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  
  // -----------
  
  SetExample(3,6);
  SampleBSpline(polyline0, nsmpl, ndegree, aKnotFlat, aCtrlPoint);
  
  // -----------
  
  glutMainLoop();
  return 0;
}


