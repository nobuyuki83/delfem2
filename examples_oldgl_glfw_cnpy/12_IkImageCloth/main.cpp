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

#include <random>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cnpy/smpl_cnpy.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/garment.h"
#include "inputs_garment.h"
#include "inputs_imgboneloc.h"
#include "delfem2/rigopt.h"
#include "delfem2/srchbv3aabb.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/caddtri_v3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/tex.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

namespace dfm2 = delfem2;

// ------------------------------------

void Draw(
    dfm2::CProjector_RigMesh& projector,
    const std::vector<dfm2::CTarget>& aTarget,
    const dfm2::opengl::CTexRGB_Rect2D& tex,
    const std::vector<dfm2::CDynTri>& aETri_Cloth,
    const std::vector<double>& aXYZ_Cloth,
    const dfm2::glfw::CViewer3& viewer)
{
  ::glEnable(GL_NORMALIZE);
  viewer.DrawBegin_oldGL();
//    Draw(aXYZ1_Body,aTri_Body,aBone,aTarget);
  tex.Draw_oldGL();
  {
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Gray(0.8f) );
    //    delfem2::opengl::DrawMeshTri3D_Edge(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
    //                                        aTri_Contact.data(), aTri_Contact.size()/3);
    // draw body
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Gray(0.8f) );
    delfem2::opengl::DrawMeshTri3D_FaceNorm(
        projector.aXYZ1_Body.data(),
        projector.aTri_Body.data(),
        projector.aTri_Body.size()/3);
    // draw cloth
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    //    delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ, aETri);
    delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ_Cloth, aETri_Cloth);
    // draw cloth
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Red() );
    delfem2::opengl::DrawMeshDynTri_FaceNorm(aETri_Cloth, aXYZ_Cloth.data());
    // draw target
    ::glEnable(GL_DEPTH_TEST);
    ::glBegin(GL_LINES);
    ::glColor3d(1,0,0);
    for(auto & it : aTarget){
      dfm2::CVec3d p = it.pos;
      dfm2::opengl::myGlVertex(p+10.0*dfm2::CVec3d(0,0,1));
      dfm2::opengl::myGlVertex(p-10.0*dfm2::CVec3d(0,0,1));
    }
    ::glEnd();
  }
  glfwSwapBuffers(viewer.window);
  glfwPollEvents();
}

int main()
{
  std::vector<dfm2::CDynTri> aETri_Cloth;
  std::vector<dfm2::CVec2d> aVec2_Cloth;
  std::vector<double> aXYZ_Cloth; // deformed vertex positions
  std::vector<unsigned int> aLine_Cloth;
  {
    std::string name_cad_in_test_input;
    double scale_adjust = 0.0;
    std::vector<unsigned int> aIESeam;
    double mesher_edge_length;
    std::vector<dfm2::CRigidTrans_2DTo3D> aRT23;
    // -------
    //    Inputs_SmplTshirt(name_cad_in_test_input,
    Inputs_SmplLtshirt(name_cad_in_test_input,
                       scale_adjust,
                       aIESeam,
                       mesher_edge_length,
                       aRT23);
    std::string path_svg = std::string(PATH_INPUT_DIR)+"/"+name_cad_in_test_input;
    std::cout << "open svg: " << path_svg << std::endl;
    delfem2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 0.001*scale_adjust);
    // -------
    dfm2::CMesher_Cad2D mesher;
    dfm2::MeshingPattern(aETri_Cloth,aVec2_Cloth,aXYZ_Cloth,aLine_Cloth,mesher,
                         aRT23,cad,aIESeam,mesher_edge_length);
  }
  std::vector<double> aXYZt_Cloth = aXYZ_Cloth;
  std::vector<double> aUVW_Cloth(aXYZ_Cloth.size(), 0.0);
  const std::vector<int> aBCFlag_Cloth(aXYZ_Cloth.size()/3, 0.0);
  std::vector<dfm2::CInfoNearest<double>> aInfoNearest_Cloth;
  // -----------------------------------
  std::vector<dfm2::CTarget> aTarget;
  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    std::string name_img_in_test_inputs;
    std::vector< std::pair<double,dfm2::CVec2d> > aBoneLoc;
    double scale = 1.0;
    BoneLocs_SmplUglysweater(name_img_in_test_inputs,
                             scale,
                             aBoneLoc);
    int width, height;
    {
      int channels;
      {
        unsigned char *img = stbi_load(
            (std::string(PATH_INPUT_DIR) + "/" + name_img_in_test_inputs).c_str(),
            &width, &height, &channels, 0);
        tex.Initialize(width, height, channels, img);
        stbi_image_free(img);
      }
      tex.max_x = -scale*width*0.5;
      tex.min_x = +scale*width*0.5;
      tex.max_y = -scale*height*0.5;
      tex.min_y = +scale*height*0.5;
      tex.z = -0.5;
    }
    aTarget.clear();
    for(auto & it : aBoneLoc){
      dfm2::CTarget t;
      t.ib = it.first;
      int iw = (int)it.second.x;
      int ih = (int)it.second.y;
      t.pos.p[0] = (double)iw/width-0.5;
      t.pos.p[1] = 0.5*height/width - (double)ih/height;
      aTarget.push_back(t);
    }
  }

  dfm2::CProjector_RigMesh projector;
  {
    std::vector<unsigned int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    std::vector<double> aW_Body;
    dfm2::cnpy::LoadSmpl_Bone(
        projector.aXYZ0_Body,
        aW_Body,
        projector.aTri_Body,
        aIndBoneParent,
        aJntRgrs,
        std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    {
      std::vector<double> aJntPos0;
      dfm2::Points3_WeighttranspPosition(
          aJntPos0,
          aJntRgrs, projector.aXYZ0_Body);
      dfm2::InitBones_JointPosition(
          projector.aBone,
          aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
    }
//    dfm2::Smpl2Rig(projector.aBone,
//        aIndBoneParent, projector.aXYZ0_Body, aJntRgrs);
    dfm2::SparsifyMatrixRow(
        projector.aSkinningSparseWeight,
        projector.aSkinningSparseIdBone,
        aW_Body.data(),
        projector.aXYZ0_Body.size()/3,
        projector.aBone.size(),
        1.0e-5);
    projector.UpdatePose(true);
  }

  std::vector< std::pair<dfm2::CVec3d,dfm2::CVec3d> > aTargetOriginPos;
  for(auto & target : aTarget){
    unsigned int ib = target.ib;
    aTargetOriginPos.emplace_back(target.pos, projector.aBone[ib].Pos() );
  }
  dfm2::CKineticDamper damper;
  dfm2::glfw::CViewer3 viewer;
  // -------------------
  // below: opengl starts
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::YTOP;
  viewer.camera.view_height = 1.0;
  dfm2::opengl::setSomeLighting();
  tex.InitGL();

  const double dt = 0.01;
  const double gravity[3] = {0.0, -0.1, 0.0};
  const double bend_stiff_ratio = 0.01;

  {
    for(int iframe=0;iframe<300;++iframe){
      dfm2::StepTime_PbdClothSim(
          aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
          aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
          projector,
          dt,gravity,bend_stiff_ratio);
      damper.Damp(aUVW_Cloth);
      Draw(projector,aTarget,tex,aETri_Cloth,aXYZ_Cloth,viewer);
      viewer.ExitIfClosed();
    }
    for(int iframe=0;iframe<1000;++iframe){
      double r = (double)(iframe)/1000.0;
      if( r > 1 ){ r = 1; }
      for(unsigned int it=0;it<aTarget.size();++it){
        aTarget[it].pos = r*aTargetOriginPos[it].first + (1-r)*aTargetOriginPos[it].second;
      }
      Solve_MinRigging(projector.aBone, aTarget);
      projector.UpdatePose(false);
      dfm2::StepTime_PbdClothSim(
          aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
          aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
          projector,
          dt,gravity,bend_stiff_ratio);
      damper.Damp(aUVW_Cloth);
      Draw(projector,aTarget,tex,aETri_Cloth,aXYZ_Cloth,viewer);
      viewer.ExitIfClosed();
    }
    for(int iframe=0;iframe<300;iframe++){
      dfm2::StepTime_PbdClothSim(
          aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
          aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
          projector,
          dt,gravity,bend_stiff_ratio);
      damper.Damp(aUVW_Cloth);
      Draw(projector,aTarget,tex,aETri_Cloth,aXYZ_Cloth,viewer);
      viewer.ExitIfClosed();
    }
  }
}
