#include <iostream>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/primitive.h"

#include "delfem2/gl2_funcs.h"
#include "delfem2/gl_color.h"

#include "../glut_funcs.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

CNav3D_GLUT window;
std::vector< std::vector<unsigned int> > aaQuad;
std::vector< std::vector<double> > aaXYZ;
const unsigned int nlevel_subdiv = 3;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  
  DrawBackground(CColor(0.2,0.7,0.7));
  
  ::glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  ::glEnable(GL_LIGHTING);
  DrawMeshQuad3D_FaceNorm(aaXYZ[nlevel_subdiv],aaQuad[nlevel_subdiv]);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  DrawMeshQuad3D_Edge(aaXYZ[nlevel_subdiv],aaQuad[nlevel_subdiv]);
  
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
    case 'a':
      break;
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
  
  aaXYZ.resize(nlevel_subdiv+1);
  aaQuad.resize(nlevel_subdiv+1);
  MeshQuad3D_CubeVox(aaXYZ[0],aaQuad[0],
                     -1,+1,  -1,+1,  -1,+1);
  
  for(unsigned int il=0;il<nlevel_subdiv;++il){
    const std::vector<double>& aXYZ0 = aaXYZ[il];
    const std::vector<unsigned int>& aQuad0 = aaQuad[il];
    std::vector<unsigned int>& aQuad1 = aaQuad[il+1];
    std::vector<int> aEdgeFace0;
    std::vector<int> psupIndQuad0, psupQuad0;
    QuadSubdiv(aQuad1,
               psupIndQuad0,psupQuad0, aEdgeFace0,
               aQuad0.data(), aQuad0.size()/4,
               aXYZ0.size()/3);
    ///////
    std::vector<double>& aXYZ1 = aaXYZ[il+1];
    SubdivisionPoints_QuadCatmullClark(aXYZ1,
                                       aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,
                                       aQuad0.data(), aQuad0.size()/4,
                                       aXYZ0.data(),  aXYZ0.size()/3);
  }
  
  window.camera.view_height = 2.0;
  
  setSomeLighting();
  
  glutMainLoop();
  return 0;
}


