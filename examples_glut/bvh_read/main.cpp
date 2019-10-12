#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/vec3.h"
#include "delfem2/vec2.h"
#include "delfem2/mshio.h"
#include "delfem2/v23m3q.h"
#include "delfem2/funcs.h"

#include "delfem2/rig_v3q.h"

#include "delfem2/gl2_funcs.h"
#include "delfem2/gl_v23.h"
#include "delfem2/gl_rig_v23q.h"

#include "../glut_cam.h"


CNav3D_GLUT nav;
std::vector<CRigBone> aBone;
std::vector<CChannel_BioVisionHierarchy> aChannelRotTransBone;
int nframe = 0;
std::vector<double> aValRotTransBone;

int iframe;
int ibone_selected = -1;
int ielem_bone_selected = -1;
bool is_animation = false;
double rad_bone_sphere = 0.1;
double rad_rot_hndlr = 1.0;

void myGlutDisplay(void)
{
  //  ::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(0.5f, 0.8f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  nav.SetGL_Camera();
  
  DrawAxis(10);
  
  DrawBone(aBone,
           ibone_selected, ielem_bone_selected,
           rad_bone_sphere, rad_rot_hndlr);
  
  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}


void myGlutIdle()
{
  if( is_animation ){
    const int nch = aChannelRotTransBone.size();
    SetPose_BioVisionHierarchy(aBone,aChannelRotTransBone,aValRotTransBone.data()+iframe*nch);
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
  nav.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  nav.glutMotion(x, y);
  if( nav.imodifier != 0 ) return;
  ////
  if( ibone_selected>=0 && ibone_selected<aBone.size() ){
    nav.SetGL_Camera();
    float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
    float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
    CVector2 sp1(nav.mouse_x, nav.mouse_y);
    CVector2 sp0(nav.mouse_x-nav.dx, nav.mouse_y-nav.dy);
    CRigBone& bone = aBone[ibone_selected];
    DragHandlerRot_Mat4(bone.rot,
                        ielem_bone_selected, sp0, sp1, bone.Mat,
                        mMV, mPj);
    UpdateBoneRotTrans(aBone);
  }
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button, state, x, y);
  /////
  nav.SetGL_Camera();
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector3 src = screenUnProjection(CVector3(nav.mouse_x,nav.mouse_y,0.0), mMV, mPj);
  CVector3 dir = screenDepthDirection(src,mMV,mPj);
  if( state == GLUT_DOWN ){
    const double wh = 1.0/mPj[5];
    PickBone(ibone_selected, ielem_bone_selected,
             aBone,
             src,dir,
             rad_rot_hndlr,
             wh*0.02);
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
  
  nav.camera.view_height = 15.0;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  
  std::string path_bvh = std::string(PATH_INPUT_DIR)+"/walk.bvh";
  std::cout << "path:" << path_bvh << std::endl;
  Read_BioVisionHierarchy(aBone,aChannelRotTransBone,nframe,aValRotTransBone,
          path_bvh);
  UpdateBoneRotTrans(aBone);
  std::cout << "nBone:" << aBone.size() << "   aCh:" << aChannelRotTransBone.size() << std::endl;
  for(unsigned int ib=0;ib<aBone.size();++ib){
    std::cout << ib << " " << aBone[ib].name << std::endl;
  }
  iframe = 0;

  
  glutMainLoop();
  return 0;
}
