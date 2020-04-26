/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief SMPL model
 * @details skinning
 */

#include <cstdlib>
#include "delfem2/rig_geo3.h"
#include "delfem2/quat.h"
#include "delfem2/mat4.h"
#include "delfem2/gizmo_geo3.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/gizmo_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/rigv3_glold.h"
#include "delfem2/opengl/tex_gl.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace dfm2 = delfem2;


int main()
{
  class CMyViewer : public delfem2::opengl::CViewer_GLFW {
  public:
    CMyViewer(){
      std::vector<int> aIndBoneParent;
      std::vector<double> aJntRgrs;
      dfm2::cnpy::LoadSmpl(aXYZ0,
                           aW,
                           aTri,
                           aIndBoneParent,
                           aJntRgrs,
                           std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
      Smpl2Rig(aBone,
               aIndBoneParent, aXYZ0, aJntRgrs);
      aXYZ1 = aXYZ0;
      gizmo_rot.size = 0.3;
    }
    void Draw(){
      ::glEnable(GL_LIGHTING);
      ::glEnable(GL_DEPTH_TEST);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
      ::glDisable(GL_DEPTH_TEST);
      delfem2::opengl::DrawBone(aBone,
                                -1, -1,
                                0.01, 1.0);
      //    dfm2::opengl::DrawJoints(aJntPos1, aIndBoneParent);
      //    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
      //    dfm2::opengl::DrawJoints(aJntPos0, aIndBoneParent);
      if( ipicked_bone != -1 ){
        assert( ipicked_bone >= 0 && ipicked_bone < aBone.size() );
        gizmo_rot.pos = aBone[ipicked_bone].Pos().Float();
        { // set quaternion
          dfm2::CMat3<double> m3;
          m3.SetMat4(aBone[ipicked_bone].affmat3Global);
          dfm2::CQuat<double> qj;
          m3.GetQuat_RotMatrix(qj.q);
          qj.CopyTo(gizmo_rot.quat);
        }
        dfm2::opengl::Draw(gizmo_rot);
      }
    }
    void SetRandomPose(){
      for(unsigned int ibone=0;ibone<aBone.size();++ibone){
        dfm2::CQuatd::Random(0.2).CopyTo(aBone[ibone].quatRelativeRot);
      }
      dfm2::CVec3d::Random().CopyToScale(aBone[0].transRelative, 0.2);
      dfm2::UpdateBoneRotTrans(aBone);
      dfm2::Skinning_LBS(aXYZ1,
                         aXYZ0, aBone, aW);
    }
    virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) {
      if( ipicked_bone != -1 ){
        assert( ipicked_bone >= 0 && ipicked_bone < aBone.size() );
        gizmo_rot.pos = aBone[ipicked_bone].Pos().Float();
        gizmo_rot.Drag(src0, src1, dir);
        {
          const int ibp = aBone[ipicked_bone].ibone_parent;
          dfm2::CMat3d m3;
          if( ibp == -1 ){ m3.SetIdentity(); }
          else{ m3.SetMat4(aBone[ibp].affmat3Global); }
          dfm2::CQuatd qp; m3.GetQuat_RotMatrix(qp.q);
          dfm2::CQuatd qg = dfm2::CQuatf(gizmo_rot.quat).Double();
          dfm2::CQuatd qj = qp.Conjugate()*qg;
          qj.CopyTo(aBone[ipicked_bone].quatRelativeRot);
        }
        dfm2::UpdateBoneRotTrans(aBone);
        dfm2::Skinning_LBS(aXYZ1,
                           aXYZ0, aBone, aW);
      }
    }
    virtual void mouse_press(const float src[3], const float dir[3])
    {
      const dfm2::CVec3d s0(src), d0(dir);
      if( ipicked_bone != -1 ){
        gizmo_rot.Pick(true, src, dir, 0.03);
        if( gizmo_rot.ielem_picked != -1 ){
          return;
        }
      }
      ipicked_bone = -1;
      for(unsigned int ib=0;ib<aBone.size();++ib){
        dfm2::CVec3d p0 = aBone[ib].Pos();
        dfm2::CVec3d p1 = dfm2::nearest_Line_Point(p0, s0, d0);
        double len01 = (p0-p1).Length();
        if( len01 < 0.03 ){
          ipicked_bone = ib;
          gizmo_rot.quat[0] = aBone[ipicked_bone].quatRelativeRot[0];
          gizmo_rot.quat[1] = aBone[ipicked_bone].quatRelativeRot[1];
          gizmo_rot.quat[2] = aBone[ipicked_bone].quatRelativeRot[2];
          gizmo_rot.quat[3] = aBone[ipicked_bone].quatRelativeRot[3];
        }
      }
    }
  public:
    std::vector<double> aXYZ0, aXYZ1;
    std::vector<double> aW;
    std::vector<unsigned int> aTri;
    std::vector<dfm2::CRigBone> aBone;
    dfm2::CGizmo_Rotation<float> gizmo_rot;
    int ipicked_bone = -1;
  } viewer;
  // -----------
  dfm2::opengl::CTexRGB_Rect2D tex;
  int width, height;
  {
    int channels;
    unsigned char *img = stbi_load((std::string(PATH_INPUT_DIR)+"/uglysweater.jpg").c_str(),
                                   &width, &height, &channels, 0);
    tex.Initialize(width, height, img, "rgb");
    double scale = 0.0020;
    delete[] img;
    tex.max_x = -scale*width*0.5;
    tex.min_x = +scale*width*0.5;
    tex.max_y = -scale*height*0.5;
    tex.min_y = +scale*height*0.5;
    tex.z = -0.5;
    std::cout << width << " " << height << std::endl;
  }
  // -------------------
  viewer.Init_oldGL();
  tex.InitGL();
  dfm2::opengl::setSomeLighting();
  while (true)
  {
    viewer.DrawBegin_oldGL();
    tex.Draw_oldGL();
    viewer.Draw();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
