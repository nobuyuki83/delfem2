/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/gizmo_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/gizmo_geo3.h"
#include "delfem2/mshio.h"
#include "delfem2/points.h"
#include "delfem2/mat4.h"
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// -------------------------------------------

int main(int argc,char* argv[])
{
  class CMyViewer : public delfem2::opengl::CViewer_GLFW {
  public:
    CMyViewer(){
      delfem2::Read_Ply(
          std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
          aXYZ,aTri);
      delfem2::Normalize_Points3(aXYZ);
      gizmo_rot.size = 0.7;
    }
    //
    void mouse_press(const float src[3], const float dir[3]) override{
      gizmo_rot.Pick(true, src, dir, 0.1);
    }
    void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) override{
      gizmo_rot.Drag(src0, src1, dir);
    }
    //
    void Draw(){
      DrawBegin_oldGL();
      {
        float r[16]; dfm2::Mat4_Quat(r, gizmo_rot.quat);
        float r0[16]; dfm2::Transpose_Mat4(r0, r);
        ::glMatrixMode(GL_MODELVIEW);
        ::glPushMatrix();
        ::glMultMatrixf(r0);
        ::glEnable(GL_LIGHTING);
        ::glColor3d(0,0,0);
        delfem2::opengl::DrawMeshTri3D_Edge(
            aXYZ.data(), aXYZ.size()/3,
            aTri.data(), aTri.size()/3);
        delfem2::opengl::DrawMeshTri3D_FaceNorm(
            aXYZ.data(),
            aTri.data(), aTri.size()/3);
        ::glMatrixMode(GL_MODELVIEW);
        ::glPopMatrix();
      }
      delfem2::opengl::Draw(gizmo_rot);
      SwapBuffers();
    }
  public:
    dfm2::CGizmo_Rotation<float> gizmo_rot;
    std::vector<double> aXYZ;
    std::vector<unsigned int> aTri;
  } viewer;
  // --------------------
  viewer.Init_oldGL();
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  // --------------------
  while(true){
    viewer.Draw();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


