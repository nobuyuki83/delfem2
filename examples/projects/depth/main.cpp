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

class CInputDepth0: public CInputDepth
{
  void Draw() const {
    ::Draw();
  }
} input;
CDepth depth;
CContextDepth depth_context;
bool is_animation = false;

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
  depth.Draw_Point(false);
  
  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}

void myGlutIdle(){
  if(is_animation){
    glBindFramebuffer(GL_FRAMEBUFFER, depth_context.id_framebuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depth_context.id_depth_render_buffer);
    depth.TakeDepthShot(input);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
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
  
  depth_context.Init(500, 500);
  
  depth.nResW = 500;
  depth.nResH = 500;
  depth.lengrid = 0.005;
  depth.depth_max = 4.0;
  depth.dirPrj = CVector3(0,0,-1);
  depth.dirWidth = CVector3(1,0,0);
  depth.orgPrj = CVector3(0,0,-2);
  depth.SetColor(1, 0, 0);
 
  glutMainLoop();
	return 0;
}


