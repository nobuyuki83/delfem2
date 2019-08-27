#include <iostream>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/emat.h"

#include "delfem2/gl_funcs.h"
#include "delfem2/gl_color.h"
#include "delfem2/glut_funcs.h"


// area of a triangle
static double TriArea2D(const double p0[], const double p1[], const double p2[]){
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

// display data
bool is_animation;

CGlutWindowManager window;

void drawObject(){ // for shadow
  ::glutSolidTeapot(1.0);
}

void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  window.SetGL_Camera();
  
  DrawBackground();
  
  
  ::drawFloorShadow(drawObject, -5, 20);
  
  ShowFPS();
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
    case 'c':
    {
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
  
  ////////////////////////
  
  window.camera.view_height = 2.0;
  setSomeLighting();
  
  {
    for(int itr=0;itr<200;++itr){
      double C[3][2];
      for(int i=0;i<6;++i){
        (&C[0][0])[i] = 10.0*(rand()/(RAND_MAX+1.0)-0.5);
      }
      double a0 = TriArea2D(C[0], C[1], C[2]);
      if( a0 < 0.1 ) continue;
      double u[3][3];
      for(int i=0;i<9;++i){
        (&u[0][0])[i] = 1.0*(rand()/(RAND_MAX+1.0)-0.5);
      }
      double thickness = (rand()+1.0)/(RAND_MAX+1.0);
      double lambda = (rand()+1.0)/(RAND_MAX+1.0);
      double myu = (rand()+1.0)/(RAND_MAX+1.0);
      double diff = Check_WdWddW_PlateBendingMITC3(C, u,
                                                   thickness,lambda,myu, 1.0e-5);
      std::cout << itr << " " << diff << std::endl;
    }
  }
  
  glutMainLoop();
  return 0;
}


