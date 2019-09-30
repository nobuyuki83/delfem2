#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#if defined(__APPLE__) && (__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/msh.h"
#include "delfem2/mshio.h"
#include "delfem2/dtri_v3.h"
#include "delfem2/dtri.h"

#include "delfem2/gl2_funcs.h"
#include "delfem2/gl_color.h"

#include "../glut_funcs.h"

void MyGlVertex3dv(CVector3& p){
  ::glVertex3d(p.x, p.y, p.z);
}

void MyGlNormal3dv(CVector3& n){
  ::glNormal3d(n.x, n.y, n.z);
}

double cur_time = 0.0;
double dt = 0.1;
std::vector<CEPo2> aPo;
std::vector<ETri> aTri;
std::vector<CVector3> aVec3;

CNav3D_GLUT nav;
bool is_animation = true;
bool is_lighting = false;

void SetNewProblem()
{
    int nnode;
    double* pXYZs = 0;
    int ntri;
    unsigned int* aTriInd = 0;
    //    Load_Ply("homer.ply" ,nnode,pXYZs, ntri,aTriInd);
    Read_Ply(std::string(PATH_INPUT_DIR)+"/arm_16k.ply" ,nnode,pXYZs, ntri,aTriInd);
    {
      double cx,cy,cz, wx,wy,wz;
      GetCenterWidth(cx,cy,cz, wx,wy,wz,
                     nnode,pXYZs);
      Translate(-cx,-cy,-cz, nnode,pXYZs);
      double wm = wx;
      wm = ( wx > wm ) ? wx : wm;
      wm = ( wy > wm ) ? wy : wm;
      wm = ( wz > wm ) ? wz : wm;
      Scale(2.0/wm,nnode,pXYZs);
    }
    aPo.resize(nnode);
    aVec3.resize(nnode);
    for(unsigned int ipo=0;ipo<aPo.size();ipo++){
      aVec3[ipo].x = pXYZs[ipo*3+0];
      aVec3[ipo].y = pXYZs[ipo*3+1];
      aVec3[ipo].z = pXYZs[ipo*3+2];
    }
    InitializeMesh(aPo, aTri,
                   aTriInd,ntri,aVec3.size());
    delete[] pXYZs;
    delete[] aTriInd;
  
  CheckTri(aTri);
  CheckTri(aPo, aTri);
  CheckTri(aPo, aTri, aVec3);
}

//////////////////////////////////////////////////////////////

void myGlutIdle(){
  if( is_animation ){
    for(unsigned int i=0;i<10;i++){
      unsigned int itri0 = (unsigned int)((rand()/(RAND_MAX+1.0))*aTri.size());
      assert( itri0 < aTri.size() );
      Collapse_ElemEdge(itri0, 0, aPo, aTri);
      CheckTri(aPo, aTri);
      if( aTri.size() <= 100 ) break;
    }
  }
  glutPostRedisplay();
}

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  ::glMatrixMode(GL_PROJECTION);
  glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  nav.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y){
  nav.glutMouse(button, state, x, y);
}

void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key, x, y);
}

void myGlutDisplay(void)
{
  ::glClearColor(0.2f, .7f, .7f,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  DrawBackground();
  nav.SetGL_Camera();
  
  GLboolean is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glEnable(GL_LIGHTING);
  GLboolean is_texture  = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_TEXTURE_2D);    
  {
    //    float gray[4] = {0.3,0.3,0.3,1};
    float gray[4] = {0.9f,0.9f,0.9f,1.f};    
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    float shine[4] = {0,0,0,0};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
    ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 127.0);
    //    ::glColor3d(1,1,1);
  }
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itri=0;itri<aTri.size();itri++){  
    const unsigned int i1 = aTri[itri].v[0];
    const unsigned int i2 = aTri[itri].v[1];
    const unsigned int i3 = aTri[itri].v[2];
    MyGlVertex3dv(aVec3[i1]);
    MyGlVertex3dv(aVec3[i2]);
    MyGlVertex3dv(aVec3[i3]);
  }
  ::glEnd();        
  
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itri=0;itri<aTri.size();itri++){  
    const unsigned int i1 = aTri[itri].v[0];
    const unsigned int i2 = aTri[itri].v[1];
    const unsigned int i3 = aTri[itri].v[2];
    MyGlVertex3dv(aVec3[i1]);     MyGlVertex3dv(aVec3[i2]);
    MyGlVertex3dv(aVec3[i2]);     MyGlVertex3dv(aVec3[i3]);
    MyGlVertex3dv(aVec3[i3]);     MyGlVertex3dv(aVec3[i1]);
  }
  ::glEnd();        
  
  
  if( is_lighting ){ ::glEnable( GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }
  if(  is_texture  ){ ::glEnable(GL_TEXTURE_2D); } 
  
  if( is_animation  ){
    cur_time += dt;
  }
  ShowFPS();
  glutSwapBuffers();
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key){
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'a':
      if( is_animation ){ is_animation = false; }
      else{ is_animation = true; }
      break;
    case ' ':
      SetNewProblem();      
      ::glMatrixMode(GL_PROJECTION);
      ::glLoadIdentity();
      //      Com::View::SetProjectionTransform(camera);
      break;
    case 'l':
      is_lighting = !is_lighting;
      if( is_lighting ){ 
        ::glEnable(GL_LIGHTING);
        std::cout << "put on light " << std::endl;
      }
      else{              
        ::glDisable(GL_LIGHTING);
      }
      break;
  }
  
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  // initialize glut
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("Surface Mesh Edge Collapse");
  
  // define call back functions
  glutIdleFunc(myGlutIdle);
  glutKeyboardFunc(myGlutKeyboard);
  glutDisplayFunc(myGlutDisplay);
  glutReshapeFunc(myGlutResize);
  glutSpecialFunc(myGlutSpecial);;
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  
  setSomeLighting();
  SetNewProblem();
  
  nav.camera.view_height = 1.5;
  
  {
    float white[3] = {1.0,1.0,1.0};
    
    ::glEnable(GL_LIGHT0);    
    float light0pos[4] = {0.5,0.5,5,0};
    ::glLightfv(GL_LIGHT0, GL_POSITION, light0pos);
    ::glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
    ////
    ::glEnable(GL_LIGHT1);    
    float light1pos[4] = {0.5,0.5,-10,0};
    ::glLightfv(GL_LIGHT1, GL_POSITION, light1pos);      
    ::glLightfv(GL_LIGHT1, GL_DIFFUSE, white);        
  }
  glutMainLoop();
  return 0;
}
