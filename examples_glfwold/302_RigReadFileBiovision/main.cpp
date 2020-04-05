/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <vector>
#include <set>
#include "delfem2/rig_v3q.h"

// ---------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/rig_v3m3q_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -------------------------------------------

int main(int argc,char* argv[])
{
  std::vector<dfm2::CRigBone> aBone;
  std::vector<dfm2::CChannel_BioVisionHierarchy> aChannelRotTransBone;
  int nframe = 0;
  std::vector<double> aValRotTransBone;

  std::string path_bvh = std::string(PATH_INPUT_DIR)+"/jump.bvh";
//  std::cout << "path:" << path_bvh << std::endl;
  Read_BioVisionHierarchy(aBone,aChannelRotTransBone,nframe,aValRotTransBone,
                          path_bvh);
  UpdateBoneRotTrans(aBone);
  std::cout << "nBone:" << aBone.size() << "   aCh:" << aChannelRotTransBone.size() << std::endl;
  for(unsigned int ib=0;ib<aBone.size();++ib){
    std::cout << ib << " " << aBone[ib].name << std::endl;
  }
  
  // -------------------------
  
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 15.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  
  // -------------------------
  
  while(!glfwWindowShouldClose(viewer.window)){
    {
      static int iframe = 0;
      const int nch = aChannelRotTransBone.size();
      SetPose_BioVisionHierarchy(aBone, aChannelRotTransBone,
                                 aValRotTransBone.data()+iframe*nch);
      iframe = (iframe+1)%nframe;
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    delfem2::opengl::DrawAxis(10);
    dfm2::opengl::DrawBone(aBone,
                           -1, -1,
                           0.1, 1.0);
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

// ----------------------------------------------------------------------

/*
 int ibone_selected = -1;
 int ielem_bone_selected = -1;
 bool is_animation = false;
 double rad_bone_sphere = 0.1;
 double rad_rot_hndlr = 1.0;
 
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
 */
