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
#include "delfem2/mat3.h"
#include "delfem2/mshio.h"
#include "delfem2/funcs.h"
#include "delfem2/quat.h"
#include "delfem2/v23m3q.h"
#include "delfem2/rig_v3q.h"

#include "delfem2/gl_funcs.h"
#include "delfem2/gl_v23.h"
#include "delfem2/gl_rig_v23q.h"
#include "delfem2/glut_funcs.h"

#include "delfem2/../../external/tinygltf/tiny_gltf.h"
#include "delfem2/../../external/io_gltf.h"


//////////////////////////////////////////////////////////////////////

CGlutWindowManager window;
std::vector<double> aXYZ0;
std::vector<unsigned int> aTri;
std::vector<double> aRigWeight;
std::vector<unsigned int> aRigJoint;
////
std::vector<double> aXYZ;
std::vector<CRigBone> aBone;
int ibone_selected = -1;
int ielem_bone_selected = -1;
double rad_bone_sphere = 0.01;
double rad_rot_hndlr = 1.0;

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
  
//  ::glDisable(GL_LIGHTING);
//  DrawMeshTri3D_Edge(aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  
  ::glEnable(GL_LIGHTING);
  DrawMeshTri3D_FaceNorm(aXYZ.data(), aTri.data(), aTri.size()/3);
  
  DrawAxis(1);
  
  
  ::glDisable(GL_DEPTH_TEST);
  DrawBone(aBone,
           ibone_selected, ielem_bone_selected,
           rad_bone_sphere, rad_rot_hndlr);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}


void myGlutIdle()
{
  if( is_animation ){
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
    CRigBone& bone = aBone[ibone_selected];
    DragHandlerRot_Mat4(bone.rot,
                        ielem_bone_selected, sp0, sp1, bone.Mat,
                        mMV, mPj);
    UpdateBoneRotTrans(aBone);
    UpdateRigSkin(aXYZ,
                  aXYZ0, aTri, aBone, aRigWeight, aRigJoint);
  }
  ////
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
             rad_rot_hndlr,
             wh*0.05);
    std::cout << ibone_selected << std::endl;
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
  
  
  {
    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;
    
    //std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.gltf";
    //bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, path_gltf);
    
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
//      std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Monster.glb";
    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/CesiumMan.glb";
  
    bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, path_gltf); // for binary glTF(.glb)
    if (!warn.empty()) { printf("Warn: %s\n", warn.c_str()); }
    if (!err.empty()) { printf("Err: %s\n", err.c_str()); }
    if (!ret) { printf("Failed to parse glTF\n"); return -1; }
    
    Print(model);
    ////
    GetMeshInfo(aXYZ0, aTri, aRigWeight, aRigJoint,
                model, 0, 0);
    if( !model.skins.empty() ){
      aBone.resize( model.skins[0].joints.size() );
      std::vector<int> mapNode2Bone( model.nodes.size(), -1);
      for(int ij=0;ij<model.skins[0].joints.size();++ij){
        int inode = model.skins[0].joints[ij];
        mapNode2Bone[inode] = ij;
      }
      SetBone(aBone,
              model, model.skins[0].skeleton, -1, mapNode2Bone);
      GetBoneBinding(aBone,
                     model);
    }
    UpdateBoneRotTrans(aBone);
    UpdateRigSkin(aXYZ,
                  aXYZ0, aTri, aBone, aRigWeight, aRigJoint);
  }
  
  window.camera.view_height = 2.0;
  window.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  setSomeLighting();
  
  glutMainLoop();
  return 0;
}
