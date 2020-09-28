/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <cmath>
#include "delfem2/gizmo_geo3.h"
#include "delfem2/mshio.h"
#include "delfem2/points.h"
// ---
#include <GLFW/glfw3.h>
#include "delfem2/opengl/gizmo_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -------------------------------------------

int main(int argc,char* argv[])
{
  class CMyViewer : public delfem2::opengl::CViewer_GLFW {
  public:
    CMyViewer(){
      delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
                        aXYZ,aTri);
      delfem2::Normalize_Points3(aXYZ);
    }
    //
    virtual void mouse_press(const float src[3], const float dir[3]){
      ga.Pick(src, dir);
    }
    virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]){
      ga.Drag(src0, src1, dir);
    }
    virtual void key_release(int key, int mods){
    }
    virtual void key_press(int key, int mods){
      if( key == GLFW_KEY_R ){ ga.igizmo_mode = 1; }
      if( key == GLFW_KEY_G ){ ga.igizmo_mode = 0; }
    }
    //
    void Draw(){
      DrawBegin_oldGL();
      delfem2::opengl::DrawAxis(1);
      {
        ::glMatrixMode(GL_MODELVIEW);
        ::glPushMatrix();
        const auto m0 = ga.Affine();
        const auto m1 = m0.Transpose();
        delfem2::opengl::MyGlMultMat(m1);
        // ------
        ::glEnable(GL_LIGHTING);
        ::glColor3d(0,0,0);
        delfem2::opengl::DrawMeshTri3D_Edge(aXYZ.data(), aXYZ.size()/3,
                                            aTri.data(), aTri.size()/3);
        delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ.data(),
                                                aTri.data(), aTri.size()/3);
        // -------
        ::glMatrixMode(GL_MODELVIEW);
        ::glPopMatrix();
      }
      delfem2::opengl::Draw(ga);
      DrawEnd_oldGL();
    }
  public:
    delfem2::CGizmo_Affine<float> ga;
    std::vector<double> aXYZ;
    std::vector<unsigned int> aTri;
  } viewer;
  // --------------------
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  // --------------------
  while(true){
    viewer.Draw();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


