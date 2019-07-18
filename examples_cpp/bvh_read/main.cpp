#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/vec3.h"
#include "delfem2/vec2.h"
#include "delfem2/mshio.h"
#include "delfem2/funcs.h"

#include "delfem2/gl_funcs.h"
#include "delfem2/v23_gl.h"
#include "delfem2/glut_funcs.h"

#include "delfem2/rigmesh.h"

CGlutWindowManager window;
std::vector<CBone_RigMsh> aBone;
std::vector<CChannel_RotTransBone_BVH> aChannelRotTransBone;
int nframe = 0;
std::vector<double> aValRotTransBone;

int iframe;
int ibone_selected = -1;
int ielem_bone_selected = -1;
bool is_animation = false;

void myGlutDisplay(void)
{
  //  ::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(0.5f, 0.8f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  window.SetGL_Camera();
  
  DrawAxis(10);
  
  DrawBone(aBone, ibone_selected, ielem_bone_selected, 0.1);
  
  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}


void myGlutIdle()
{
  if( is_animation ){
    const int nch = aChannelRotTransBone.size();
    SetRotTransBVH(aBone,aChannelRotTransBone,aValRotTransBone.data()+iframe*nch);
    iframe = (iframe+1)%nframe;
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
  if( window.imodifier != 0 ) return;
  ////
  if( ibone_selected>=0 && ibone_selected<aBone.size() ){
    window.SetGL_Camera();
    float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
    float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
    CVector2 sp1(window.mouse_x, window.mouse_y);
    CVector2 sp0(window.mouse_x-window.dx, window.mouse_y-window.dy);
    CBone_RigMsh& bone = aBone[ibone_selected];
    DragHandlerRot(bone.quat_joint,
                   ielem_bone_selected, sp0, sp1, bone.pos,
                   mMV, mPj);
    UpdateBoneRotTrans(aBone);
  }
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  window.glutMouse(button, state, x, y);
  /////
  window.SetGL_Camera();
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector3 src = screenUnProjection(CVector3(window.mouse_x,window.mouse_y,0.0), mMV, mPj);
  CVector3 dir = screenDepthDirection(src,mMV,mPj);
  if( state == GLUT_DOWN ){
    const double wh = 1.0/mPj[5];
    std::cout << wh << std::endl;
    PickBone(ibone_selected, ielem_bone_selected,
             aBone,
             src,dir,
             15,
             wh*0.04);
  }
  else{
  }
  /////
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
    {
      is_animation = !is_animation;
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
  
  window.camera.view_height = 15.0;
  window.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  
  std::string path_bvh = "../test_inputs/walk.bvh";
  std::cout << "path:" << path_bvh << std::endl;
  ReadBVH(aBone,aChannelRotTransBone,nframe,aValRotTransBone,
          path_bvh);
  std::cout << "nBone:" << aBone.size() << "   aCh:" << aChannelRotTransBone.size() << std::endl;
  for(int ib=0;ib<aBone.size();++ib){
    std::cout << ib << " " << aBone[ib].name << std::endl;
  }
  iframe = 0;

  
  glutMainLoop();
  return 0;
}
