/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief edit rigged shape
 * @details rotation for each bone and translation for root bone
 */

#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/mat4.h"
#include <cstdlib>
#include <fstream>

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/rigv3_glold.h"
#include "delfem2/opengl/tex_gl.h"
#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace dfm2 = delfem2;

// ---------------------------

int main()
{
  class CMyViewer : public delfem2::opengl::CViewer_GLFW {
  public:
    CMyViewer(){
      std::vector<double> aW;
      std::vector<int> aIndBoneParent;
      std::vector<double> aJntRgrs;
      dfm2::cnpy::LoadSmpl_Bone(
          aXYZ0,
          aW,
          aTri,
          aIndBoneParent,
          aJntRgrs,
          std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
      {
        std::vector<double> aJntPos0;
        dfm2::Points3_WeighttranspPosition(
            aJntPos0,
            aJntRgrs, aXYZ0);
        dfm2::InitBones_JointPosition(
            aBone,
            aIndBoneParent, aJntPos0);
      }
      dfm2::SparsifyMatrixRow(
          aWeightRigSparse, aIdBoneRigSparse,
          aW.data(), aXYZ0.size()/3, aBone.size(),
          1.0e-5);
      aXYZ1 = aXYZ0;
      gizmo.SetSize(0.3);
    }
    void Draw(){
      ::glEnable(GL_LIGHTING);
      ::glEnable(GL_DEPTH_TEST);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
      ::glDisable(GL_DEPTH_TEST);
      delfem2::opengl::DrawBone(aBone,
                                -1, -1,
                                0.01, 1.0);
      dfm2::opengl::Draw(gizmo,aBone);

    }
    void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) override {
      bool is_edited = gizmo.Drag(aBone,
                                  src0, src1, dir);
      if( is_edited ){
        dfm2::SkinningSparse_LBS(aXYZ1,
            aXYZ0, aBone, aWeightRigSparse,aIdBoneRigSparse);
      }
    }
    void mouse_press(const float src[3], const float dir[3]) override {
      gizmo.Pick(src, dir, aBone);
    }
    void key_press(int key, int mods) override{
      if( key == GLFW_KEY_G ){ gizmo.SetMode(dfm2::CGizmo_Rig<float>::MODE_EDIT::TRNSL); }
      if( key == GLFW_KEY_R ){ gizmo.SetMode(dfm2::CGizmo_Rig<float>::MODE_EDIT::ROT); }
      if( key == GLFW_KEY_S ){
        std::ofstream fout("pose.txt");
        for(const auto & bone: aBone){
          const double* q = bone.quatRelativeRot;
          fout << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << std::endl;
        }
        {
          const double* t = aBone[0].transRelative;
          fout << t[0] << " " << t[1] << " " << t[2] << std::endl;
        }
      }
    }
  public:
    std::vector<double> aXYZ0, aXYZ1;
    std::vector<unsigned int> aTri;
    std::vector<double> aWeightRigSparse;
    std::vector<unsigned int> aIdBoneRigSparse;
    std::vector<dfm2::CRigBone> aBone;
    dfm2::CGizmo_Rig<float> gizmo;
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
  }
}
