#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include "delfem2/vec3.h"
#include "delfem2/vec2.h"
#include "delfem2/mat3.h"
#include "delfem2/mshio.h"
#include "delfem2/funcs.h"
#include "delfem2/quat.h"
//
#include "delfem2/v23m3q.h"
#include "delfem2/rig_v3q.h"
//
#include "delfem2/../../external/io_gltf.h"

// ------------------

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_v23.h"
#include "delfem2/opengl/gl_rig_v23q.h"

// ---------------------------------

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<double> aRigWeight;
  std::vector<unsigned int> aRigJoint;
  std::vector<double> aXYZ;
  std::vector<CRigBone> aBone;
  {
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
//      std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Monster.glb";
    std::string path_glb = std::string(PATH_INPUT_DIR)+"/CesiumMan.glb";
    CGLTF gltf;
    gltf.Read(path_glb);
    gltf.Print();
    gltf.GetMeshInfo(aXYZ0, aTri, aRigWeight, aRigJoint, 0,0);
    gltf.GetBone(aBone, 0);
  }
  
  UpdateBoneRotTrans(aBone);
  aXYZ = aXYZ0;
  UpdateRigSkin(aXYZ.data(),
                aXYZ0.data(), aXYZ0.size()/3,
                aTri.data(), aTri.size()/3,
                aBone, aRigWeight.data(), aRigJoint.data());

  // --------------
  // opengl starts here
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_GLold();
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  
  while(!glfwWindowShouldClose(viewer.window)){
    // --------------------
    viewer.DrawBegin_Glold();
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ.data(), aTri.data(), aTri.size()/3);
    delfem2::opengl::DrawAxis(1);    
    
    ::glDisable(GL_DEPTH_TEST);
    DrawBone(aBone,
             -1, -1,
             0.01, 1.0);
    ::glEnable(GL_DEPTH_TEST);
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}






/*
 int ibone_selected = -1;
 int ielem_bone_selected = -1;
 double rad_bone_sphere = 0.01;
 double rad_rot_hndlr = 1.0;
 
 bool is_animation = false;
 
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
 UpdateRigSkin(aXYZ.data(),
 aXYZ0.data(), aXYZ.size()/3,
 aTri.data(), aTri.size()/3,
 aBone, aRigWeight.data(), aRigJoint.data());
 }
 ////
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
 */

