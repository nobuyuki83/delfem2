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
#include "delfem2/msh.h"

#include "delfem2/color_gl.h"
#include "delfem2/funcs_glut.h"

//////////////////////////////////////////////////////////////////////////////////////

std::vector<CVector3> aXYZ;
std::vector<int> aTri;

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
    for(int ixyz=0;ixyz<(int)aXYZ.size();ixyz++){
      myGlVertex3d(aXYZ[ixyz]);
    }
    ::glEnd();
    
    const int nTri = (int)aTri.size()/3;
    /////
    ::glColor4d(1,1,1,0.8);
    ::glBegin(GL_TRIANGLES);
    for(int itri=0;itri<(int)nTri;itri++){
      const int i1 = aTri[itri*3+0];
      const int i2 = aTri[itri*3+1];
      const int i3 = aTri[itri*3+2];
      myGlVertex3d(i1,aXYZ);
      myGlVertex3d(i2,aXYZ);
      myGlVertex3d(i3,aXYZ);
    }
    ::glEnd();
    ::glColor3d(0,0,0);
    ::glBegin(GL_LINES);
    for(int itri=0;itri<(int)nTri;itri++){
      const unsigned int i1 = aTri[itri*3+0];
      const unsigned int i2 = aTri[itri*3+1];
      const unsigned int i3 = aTri[itri*3+2];
      myGlVertex3d(i1,aXYZ);
      myGlVertex3d(i2,aXYZ);
      myGlVertex3d(i2,aXYZ);
      myGlVertex3d(i3,aXYZ);
      myGlVertex3d(i3,aXYZ);
      myGlVertex3d(i1,aXYZ);
    }
    ::glEnd();
  }
  
  ShowFPS();
  
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
      int nXYZ = 100;
      aXYZ.resize(nXYZ);
      for(int ixyz=0;ixyz<nXYZ;ixyz++){
        aXYZ[ixyz].x = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
        aXYZ[ixyz].y = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
        aXYZ[ixyz].z = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
      }
      ConvexHull(aTri,aXYZ);
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
  
  win.camera.view_height = 1.5;
  
  int nXYZ = 100;
  aXYZ.resize(nXYZ);
  for(int ixyz=0;ixyz<nXYZ;ixyz++){
    aXYZ[ixyz].x = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
    aXYZ[ixyz].y = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
    aXYZ[ixyz].z = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
  }
  
  ConvexHull(aTri,aXYZ);
  
  glutMainLoop();
  return 0;
}
