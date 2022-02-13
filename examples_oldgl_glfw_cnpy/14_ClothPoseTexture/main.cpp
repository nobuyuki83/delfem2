/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief putting cloth on SMPL modl
 * @details
 */

#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "inputs_garment.h"
#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/fem_quadratic_bending.h"
#include "delfem2/garment.h"
#include "delfem2/pbd_geo3.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/srch_bv3_aabb.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_points.h"
#include "delfem2/cad2_io_svg.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/caddtri_v3.h"
#include "delfem2/opengl/tex.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

namespace dfm2 = delfem2;

void Draw(
    const std::vector<dfm2::CDynTri>& aETri_Cloth,
    const std::vector<double>& aXYZ_Cloth,
    const std::vector<dfm2::CVec2d>& aVec2,
    const dfm2::CProjector_RigMesh& projector_smpl,
    const dfm2::opengl::CTexRGB& tex,
    dfm2::glfw::CViewer3& viewer)
{
  viewer.DrawBegin_oldGL();

  ::glDisable(GL_TEXTURE_2D);
  ::glEnable(GL_NORMALIZE);
  ::glEnable(GL_LIGHTING);
  dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Gray(0.8f) );
  //    delfem2::opengl::DrawMeshTri3D_Edge(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
  //                                        aTri_Contact.data(), aTri_Contact.size()/3);
  // draw body
  ::glEnable(GL_LIGHTING);
  dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Gray(0.8f) );
  delfem2::opengl::DrawMeshTri3D_FaceNorm(
      projector_smpl.aXYZ1_Body.data(),
      projector_smpl.aTri_Body.data(),
      projector_smpl.aTri_Body.size()/3);

  // draw cloth edge
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  //    delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ, aETri);
  delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ_Cloth, aETri_Cloth);

  // draw cloth face
  ::glEnable(GL_TEXTURE_2D);
  ::glBindTexture(GL_TEXTURE_2D, tex.id_tex);

  ::glEnable(GL_LIGHTING);
  dfm2::opengl::myGlColorDiffuse( dfm2::CColor::White() );
//  delfem2::opengl::DrawMeshDynTri_FaceNorm(aETri_Cloth, aXYZ_Cloth.data());
  delfem2::opengl::DrawMeshDynTri_FaceNormTex(aETri_Cloth, aXYZ_Cloth.data(), aVec2);

  glfwSwapBuffers(viewer.window);
  glfwPollEvents();
  viewer.ExitIfClosed();
}

// --------------------
int main()
{
  const double dt = 0.01;
  const double gravity[3] = {0.0, -0.1, 0.0};
  const double mesher_edge_length = 0.02;
  const double bend_stiffness_ratio = 0.01;

  // -----------------------------

  dfm2::opengl::CTexRGB tex;
  int width, height;
  { // texture imput data
    //double scale = 0.01;
    std::string name_img_in_test_inputs = "lenna.png";
    // -----
    int channels;
    unsigned char *img = stbi_load(
        (std::string(PATH_INPUT_DIR)+"/"+name_img_in_test_inputs).c_str(),
        &width, &height, &channels, 0);
    tex.Initialize(width, height, channels, img);
    stbi_image_free(img);
  }

  std::vector<dfm2::CDynTri> aETri_Cloth;
  std::vector<dfm2::CVec2d> aVec2_Cloth;
  std::vector<double> aXYZ_Cloth; // deformed vertex positions
  std::vector<unsigned int> aLine_Cloth;
  { // prepare clothing input data
    dfm2::CMesher_Cad2D mesher;
    std::string name_cad_in_test_input;
    double scale_adjust = 0.0;
    std::vector<unsigned int> aIESeam;
    double mesher_edge_length0;
    std::vector<dfm2::CRigidTrans_2DTo3D> aRT23;
    // -------
    dfm2::Inputs_SmplTshirt2(
        name_cad_in_test_input,
        scale_adjust,
        aIESeam,
        mesher_edge_length0,
        aRT23);
    std::string path_svg = std::string(PATH_INPUT_DIR)+"/"+name_cad_in_test_input;
    std::cout << "open svg: " << path_svg << std::endl;
    delfem2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 0.001*scale_adjust);
    // -------
    dfm2::MeshingPattern(
        aETri_Cloth,aVec2_Cloth,aXYZ_Cloth,aLine_Cloth,mesher,
        aRT23,cad,aIESeam,mesher_edge_length);
  }
  std::vector<double> aXYZt_Cloth = aXYZ_Cloth;
  std::vector<double> aUVW_Cloth(aXYZ_Cloth.size(), 0.0);
  const std::vector<int> aBCFlag_Cloth(aXYZ_Cloth.size()/3, 0);

  // ----------
  dfm2::CProjector_RigMesh body_smpl;
  {
    std::vector<double> aW_Body;
    std::vector<unsigned int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl_Bone(
        body_smpl.aXYZ0_Body,
        aW_Body,
        body_smpl.aTri_Body,
        aIndBoneParent,
        aJntRgrs,
        std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    {
      std::vector<double> aJntPos0;
      dfm2::Points3_WeighttranspPosition(
          aJntPos0,
          aJntRgrs, body_smpl.aXYZ0_Body);
      dfm2::InitBones_JointPosition(
          body_smpl.aBone,
          aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
    }
    dfm2::SparsifyMatrixRow(
        body_smpl.aSkinningSparseWeight,
        body_smpl.aSkinningSparseIdBone,
        aW_Body.data(),
        body_smpl.aXYZ0_Body.size() / 3,
        body_smpl.aBone.size(),
        1.0e-5);
  }

  body_smpl.UpdatePose(true);

  std::vector< dfm2::CQuatd > aQuatTarget;
  {
    std::ifstream fin(std::string(PATH_INPUT_DIR)+"/pose_smpl1.txt");
    for(unsigned int ib=0; ib < body_smpl.aBone.size(); ++ib){
      double d0,d1,d2,d3;
      fin >> d0 >> d1 >> d2 >> d3;
      auto q = dfm2::CQuatd(d0,d1,d2,d3);
      q.SetSmallerRotation();
      aQuatTarget.push_back( q );
    }
    {
      dfm2::CVec3d posTarget;
      double d0,d1,d2;
      fin >> d0 >> d1 >> d2;
      posTarget = dfm2::CVec3d(d0,d1,d2);
    }
  }
  dfm2::CKineticDamper damper;
  dfm2::glfw::CViewer3 viewer;
  // -----------
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  // viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::YTOP;
  dfm2::opengl::setSomeLighting();
  tex.InitGL();

  for(int iframe = 0; iframe< 300; ++iframe){
    dfm2::StepTime_PbdClothSim(
        aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
        aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
        body_smpl,
        dt, gravity, bend_stiffness_ratio);
    damper.Damp(aUVW_Cloth);
    Draw(aETri_Cloth, aXYZ_Cloth, aVec2_Cloth, body_smpl, tex, viewer);
  }
  for(int iframe = 0; iframe<1200; ++iframe){
    double r = (double)(iframe)/1200;
    if( r > 1 ){ r = 1; }
    for(unsigned int ib=0; ib < body_smpl.aBone.size(); ++ib){
      dfm2::CQuatd q = dfm2::SphericalLinearInterp( dfm2::CQuatd::Identity(), aQuatTarget[ib], r);
      q.normalize();
      q.CopyTo(body_smpl.aBone[ib].quatRelativeRot);
    }
    body_smpl.UpdatePose(iframe % 100 == 0);
    dfm2::StepTime_PbdClothSim(
        aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
        aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
        body_smpl,
        dt, gravity, bend_stiffness_ratio);
    damper.Damp(aUVW_Cloth);
    Draw(aETri_Cloth, aXYZ_Cloth, aVec2_Cloth, body_smpl, tex, viewer);
  }
  for(int iframe=0;iframe<300;++iframe){
    dfm2::StepTime_PbdClothSim(
        aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
        aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
        body_smpl,
        dt, gravity, bend_stiffness_ratio);
    damper.Damp(aUVW_Cloth);
    Draw(aETri_Cloth, aXYZ_Cloth, aVec2_Cloth, body_smpl, tex, viewer);
  }

  // output
  dfm2::Write_Obj_UniformMesh(
      std::string(PATH_OUTPUT_DIR) + "/body.obj",
      body_smpl.aXYZ1_Body.data(), body_smpl.aXYZ1_Body.size() / 3,
      body_smpl.aTri_Body.data(), body_smpl.aTri_Body.size() / 3, 3);

  {
    const unsigned int np = aXYZ_Cloth.size() / 3;
    const unsigned int nt = aETri_Cloth.size();
    std::ofstream fout(std::string(PATH_OUTPUT_DIR) + "/cloth.obj", std::ofstream::out);
    for (unsigned int ip = 0; ip < np; ip++) {
      fout << "v ";
      fout << aXYZ_Cloth[ip * 3 + 0] << " ";
      fout << aXYZ_Cloth[ip * 3 + 1] << " ";
      fout << aXYZ_Cloth[ip * 3 + 2] << std::endl;
    }
    for (unsigned int ip = 0; ip < np; ip++) {
      fout << "vt ";
      fout << aVec2_Cloth[ip].x << " ";
      fout << aVec2_Cloth[ip].y << " ";
      fout << 0.0 << std::endl;
    }
    for (unsigned int itri = 0; itri < nt; itri++) {
      const unsigned int i0 = aETri_Cloth[itri].v[0] + 1;
      const unsigned int i1 = aETri_Cloth[itri].v[1] + 1;
      const unsigned int i2 = aETri_Cloth[itri].v[2] + 1;
      fout << "f ";
      fout << i0 << "/" << i0 << " ";
      fout << i1 << "/" << i1 << " ";
      fout << i2 << "/" << i2 << std::endl;
    }
  }

}
