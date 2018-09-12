#include <iostream>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include "delfem2/funcs_glut.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/funcs_glew.h"
#include "delfem2/v23_gl.h"
#include "delfem2/color_gl.h"

#include "delfem2/depth_v3_gl.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

double cur_time = 0.0;

void Draw(){
  ::glRotated(+cur_time, 1,0,0);
  glutSolidTeapot(1.0);
  ::glRotated(-cur_time, 1,0,0);
}

CGPUSampler sampler;
CFrameBufferManager fbm;
bool is_animation = true;

bool is_depth = false;
CGlutWindowManager window;

void myGlutDisplay(void)
{
 	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  window.SetGL_Camera();
  
  DrawBackground( CColor(0.2,0.7,0.7) );
//  ::glDisable(GL_LIGHTING);
  ::glEnable(GL_LIGHTING);
  
  ::glColor3d(1,1,1);
  Draw();
  
  ///////

  glPointSize(3);
  sampler.Draw();
  
  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}

void myGlutIdle(){
  if(is_animation){
    fbm.Start();
    sampler.Start();
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1,1,1);
    ::glEnable(GL_LIGHTING);
    Draw();
    sampler.End();
    sampler.LoadTex(); // move the sampled image to a texture
    fbm.End();
    cur_time += 1;
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
  window.glutSpecial(Key, x, y);
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  window.glutMotion(x, y);
	::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  window.glutMouse(button, state, x, y);
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
      is_depth = !is_depth;
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
	glutInitWindowSize(300, 300);
 	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("3D View");
	glutDisplayFunc(myGlutDisplay);
	glutIdleFunc(myGlutIdle);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
  
  ////////////////////////
  
  window.camera.view_height = 2.0;
  
  setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);
  
  fbm.Init(512, 512, true,true);
  
  int nres = 128;
  double elen = 0.02;
  sampler.Init(nres, nres, true,true);
  sampler.SetCoord(elen, 4.0,
                   CVector3(-nres*elen*0.5,nres*elen*0.5,-2).stlvec(),
                   CVector3(0,0,-1).stlvec(),
                   CVector3(1,0,0).stlvec() );
  sampler.SetColor(1, 0, 0);
  sampler.draw_len_axis = 1.0;
 
  glutMainLoop();
	return 0;
}


