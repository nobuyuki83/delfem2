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
#include "delfem2/garment.h"
#include "delfem2/srchuni_v3.h"
#include "delfem2/pbd_geo3.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/srchbv3aabb.h"
#include "delfem2/mshmisc.h"
#include "delfem2/cad2_io_svg.h"
#include "delfem2/points.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/caddtri_v3.h"

namespace dfm2 = delfem2;

void Draw(
    const std::vector<dfm2::CDynTri>& aETri_Cloth,
    const std::vector<double>& aXYZ_Cloth,
    const dfm2::CProjector_RigMesh& projector_smpl,
    dfm2::glfw::CViewer3& viewer)
{
  viewer.DrawBegin_oldGL();
  {
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
    // draw cloth
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    //    delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ, aETri);
    delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ_Cloth, aETri_Cloth);
    // draw cloth
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Red() );
    delfem2::opengl::DrawMeshDynTri_FaceNorm(aETri_Cloth, aXYZ_Cloth.data());
  }
  glfwSwapBuffers(viewer.window);
  glfwPollEvents();
}

// --------------------
int main()
{
  const double dt = 0.01;
  const double gravity[3] = {0.0, -0.1, 0.0};
  const double mesher_edge_length = 0.02;
  const double bend_stiffness_ratio = 0.01;

  // -----------------------------
  // below: input data
  std::vector<dfm2::CDynTri> aETri_Cloth;
  std::vector<dfm2::CVec2d> aVec2_Cloth;
  std::vector<double> aXYZ0_Cloth; // deformed vertex positions
  std::vector<unsigned int> aLine_Cloth;
  {
    dfm2::CMesher_Cad2D mesher;
    std::string name_cad_in_test_input;
    double scale_adjust = 0.0;
    std::vector<unsigned int> aIESeam;
    double mesher_edge_length0;
    std::vector<dfm2::CRigidTrans_2DTo3D> aRT23;
    // -------
    //    Inputs_SmplTshirt(name_cad_in_test_input,
    Inputs_SmplLtshirt(
        name_cad_in_test_input,
        scale_adjust,
        aIESeam,
        mesher_edge_length0,
        aRT23);
    std::string path_svg = std::string(PATH_INPUT_DIR)+"/"+name_cad_in_test_input;
    std::cout << "open svg: " << path_svg << std::endl;
    delfem2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(
        cad, path_svg, 0.001*scale_adjust);
    // -------
    dfm2::MeshingPattern(
        aETri_Cloth,aVec2_Cloth,aXYZ0_Cloth,aLine_Cloth,mesher,
        aRT23,cad,aIESeam,mesher_edge_length);
  }

  // ----------
  dfm2::CProjector_RigMesh body;
  {
    std::vector<double> aW_Body;
    std::vector<unsigned int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl_Bone(
        body.aXYZ0_Body,
        aW_Body,
        body.aTri_Body,
        aIndBoneParent,
        aJntRgrs,
        std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    {
      std::vector<double> aJntPos0;
      dfm2::Points3_WeighttranspPosition(
          aJntPos0,
          aJntRgrs, body.aXYZ0_Body);
      dfm2::InitBones_JointPosition(
          body.aBone,
          aIndBoneParent.size(), aIndBoneParent.data(), aJntPos0.data());
    }
    dfm2::SparsifyMatrixRow(
        body.aSkinningSparseWeight,
        body.aSkinningSparseIdBone,
        aW_Body.data(),
        body.aXYZ0_Body.size() / 3,
        body.aBone.size(),
        1.0e-5);
  }

  std::vector< dfm2::CQuatd > aQuatTarget;
  {
    std::ifstream fin(std::string(PATH_INPUT_DIR)+"/pose_smpl1.txt");
    for(unsigned int ib=0; ib < body.aBone.size(); ++ib){
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
  // -------------------------------
  // below: opengl related codes
  dfm2::glfw::CViewer3 viewer;
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  //viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::YTOP;
  dfm2::opengl::setSomeLighting();

  {
    { // pose initializeation
      for(auto & ib : body.aBone){
        dfm2::CQuatd::Identity().CopyTo(ib.quatRelativeRot);
      }
      body.UpdatePose(true);
    }
    std::vector<double> aXYZ_Cloth = aXYZ0_Cloth;
    std::vector<double> aXYZt_Cloth = aXYZ0_Cloth;
    std::vector<double> aUVW_Cloth(aXYZ0_Cloth.size(), 0.0);
    const std::vector<int> aBCFlag_Cloth(aXYZ0_Cloth.size()/3, 0);
    dfm2::CKineticDamper damper;
    for(int iframe = 0; iframe< 100; ++iframe){
      dfm2::StepTime_PbdClothSim(
          aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
          aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
          body,
          dt, gravity, bend_stiffness_ratio);
      damper.Damp(aUVW_Cloth);
      Draw(aETri_Cloth, aXYZ_Cloth, body, viewer);
      viewer.ExitIfClosed();
    }
    std::vector<double> aSkinningSparseWeight_Cloth;
    std::vector<unsigned int> aSkinningSparseIdBone_Cloth;
    { // rigging the clothing mesh using the body's skeleton
      const unsigned int nb = body.aBone.size();
      const unsigned int npb = body.aXYZ0_Body.size() / 3;
      const unsigned int npc = aXYZ_Cloth.size()/3;
      std::vector<double> aW_Cloth;
      aW_Cloth.assign(npc*nb, 0.0);
      for(unsigned int ipc=0;ipc<npc;++ipc){
        dfm2::PointOnSurfaceMesh<double> pesb = dfm2::Nearest_Point_MeshTri3D(
            dfm2::CVec3d(aXYZ_Cloth.data()+ipc*3),
            body.aXYZ1_Body,
            body.aTri_Body );
        const double aRb[3] = { pesb.r0, pesb.r1, 1-pesb.r0-pesb.r1 };
        const unsigned int nbone_rigbody = body.aSkinningSparseWeight.size() / npb;
        for(unsigned int inob=0;inob<3;++inob){
          const unsigned int ipb0 = body.aTri_Body[pesb.itri * 3 + inob];
          for(unsigned int iib=0;iib<nbone_rigbody;++iib){
            const unsigned int ib0 = body.aSkinningSparseIdBone[ipb0 * nbone_rigbody + iib];
            const double wb0 = body.aSkinningSparseWeight[ipb0 * nbone_rigbody + iib];
            aW_Cloth[ipc*nb+ib0] += wb0*aRb[inob];
          }
        }
      }
      dfm2::SparsifyMatrixRow(
          aSkinningSparseWeight_Cloth,
          aSkinningSparseIdBone_Cloth,
          aW_Cloth.data(),
          aXYZ_Cloth.size() / 3,
          body.aBone.size(),
          1.0e-2);
      std::cout << aSkinningSparseIdBone_Cloth.size() / (aXYZ_Cloth.size() / 3) << std::endl;
    }
    const std::vector<double> aXYZ1_Cloth = aXYZ_Cloth;
    for(int iframe = 0; iframe<30; ++iframe){
      double r = (double)(iframe)/30;
      if( r > 1 ){ r = 1; }
      for(unsigned int ib=0; ib < body.aBone.size(); ++ib){
        dfm2::CQuatd q = dfm2::SphericalLinearInterp( dfm2::CQuatd::Identity(), aQuatTarget[ib], r);
        q.normalize();
        q.CopyTo(body.aBone[ib].quatRelativeRot);
      }
      body.UpdatePose(iframe % 100 == 0);
      SkinningSparse_LBS(aXYZ_Cloth,
                        aXYZ1_Cloth, body.aBone, aSkinningSparseWeight_Cloth, aSkinningSparseIdBone_Cloth);
      dfm2::StepTime_PbdClothSim(
          aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
          aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
          body,
          dt, gravity, bend_stiffness_ratio);
      damper.Damp(aUVW_Cloth);
      Draw(aETri_Cloth, aXYZ_Cloth, body, viewer);
      viewer.ExitIfClosed();
    }
    for(int iframe=0;iframe<100;++iframe){
      dfm2::StepTime_PbdClothSim(
          aXYZ_Cloth, aXYZt_Cloth, aUVW_Cloth,
          aBCFlag_Cloth, aETri_Cloth, aVec2_Cloth, aLine_Cloth,
          body,
          dt, gravity, bend_stiffness_ratio);
      damper.Damp(aUVW_Cloth);
      Draw(aETri_Cloth, aXYZ_Cloth, body, viewer);
      viewer.ExitIfClosed();
    }
  }
}
