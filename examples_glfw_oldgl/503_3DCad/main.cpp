#include <iostream>
#include <sstream>
#include <math.h>
#include <complex>
#include <set>
#include <stack>
#include "delfem2/mat3.h"
#include "delfem2/cad3d.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

// ----------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_v23.h"

// ------------------------------

int main(int argc,char* argv[])
{
  class CCAD3DViewer : public CViewer_GLFW {
  public:
    CCAD3DViewer(){
      cad.Initialize_Sphere();
      nav.camera.view_height = 1.5;
      nav.camera.camera_rot_mode = CAMERA_ROT_YTOP;
      cad.imode_edit = CCad3D::EDIT_MOVE;
    }
    void Draw(){
      this->DrawBegin_Glold();
      cad.DrawFace_RightSelected(false);
      cad.DrawVtxEdgeHandler(nav.camera.view_height);
      this->DrawEnd_oldGL();
    }
    virtual void mouse_press(const float src[3], const float dir[3]){
      const CVector3 src_pick(src), dir_pick(dir);
      float mMV[16], mPrj[16]; nav.Matrix_MVP(mMV, mPrj, this->window);
      cad.MouseDown(src_pick, dir_pick,
                    CVector2(nav.mouse_x,nav.mouse_y),
                    mMV,mPrj,
                    nav.camera.view_height);
    }
    virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]){
      CVector2 sp0(nav.mouse_x-nav.dx, nav.mouse_y-nav.dy);
      CVector2 sp1(nav.mouse_x, nav.mouse_y);
      const CVector3 src_pick(src1), dir_pick(dir);
      float mMV[16], mPrj[16]; nav.Matrix_MVP(mMV, mPrj, this->window);
      cad.MouseMotion(src_pick,dir_pick, sp0,sp1, mMV, mPrj);
    }
  public:
    CCad3D cad;
  };
  // -------------
  CCAD3DViewer viewer;
  viewer.Init_GLold();
  delfem2::opengl::setSomeLighting();
  while(!glfwWindowShouldClose(viewer.window)){
    viewer.Draw();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


