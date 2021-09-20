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
#include "delfem2/pbd_geo3.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/srchbv3aabb.h"
#include "delfem2/mshmisc.h"
#include "delfem2/cad2_io_svg.h"
#include "delfem2/openglstb/glyph.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/opengl/old/cad2dtriv2.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/caddtri_v3.h"

namespace dfm2 = delfem2;

void Draw(
    const std::vector<dfm2::CDynTri>& aETri_Cloth,
    std::vector<unsigned int> aLine_Cloth,
    const std::vector<double>& aXYZ_Cloth,
    [[maybe_unused]] const std::vector<dfm2::CVec2d>& aVec2,
    const dfm2::CProjector_RigMesh& projector_smpl,
    dfm2::glfw::CViewer3& viewer)
{
  viewer.DrawBegin_oldGL();

  { // draw seam
    ::glDisable(GL_TEXTURE_2D);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,1,1);
    const unsigned int nline = aLine_Cloth.size()/2;
    ::glBegin(GL_LINES);
    for(unsigned int il=0;il<nline;++il) {
      const unsigned int ip0 = aLine_Cloth[il * 2 + 0];
      const unsigned int ip1 = aLine_Cloth[il * 2 + 1];
      ::glVertex3dv(aXYZ_Cloth.data()+ip0*3);
      ::glVertex3dv(aXYZ_Cloth.data()+ip1*3);
    }
    ::glEnd();
  }

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

  ::glEnable(GL_LIGHTING);
  dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Red() );
  delfem2::opengl::DrawMeshDynTri_FaceNorm(aETri_Cloth, aXYZ_Cloth.data());
//  delfem2::opengl::DrawMeshDynTri_FaceNormTex(aETri_Cloth, aXYZ_Cloth.data(), aVec2);

  glfwSwapBuffers(viewer.window);
}

// --------------------
int main()
{
  const double mesher_edge_length = 0.02;

  // -----------------------------

  std::vector<dfm2::CDynTri> aETri_Cloth;
  std::vector<dfm2::CVec2d> aVec2_Cloth;
  std::vector<double> aXYZ_Cloth; // deformed vertex positions
  std::vector<unsigned int> aLine_Cloth;
  delfem2::CCad2D cad;
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

  delfem2::openglstb::CGlyph glyph(std::string(PATH_INPUT_DIR)+"/myFont.png");
  glyph.ParseGlyphInfo(std::string(PATH_INPUT_DIR)+"/myFont.fnt");

   // -----------
  dfm2::glfw::CViewer3 viewer3;
  dfm2::glfw::CViewer2 viewer2;
  dfm2::glfw::InitGLOld();
  viewer3.InitGL();
  // viewer3.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::YTOP;
  viewer3.projection.view_height = 1.0;
  dfm2::opengl::setSomeLighting();
  viewer2.InitGL();
  ::glfwMakeContextCurrent(viewer2.window);
  glyph.InitGL();

  for(;;) {
    ::glfwMakeContextCurrent(viewer3.window);
    Draw(aETri_Cloth, aLine_Cloth, aXYZ_Cloth, aVec2_Cloth, body_smpl, viewer3);

    viewer2.DrawBegin_oldGL();
    delfem2::opengl::Draw_CCad2D(cad);
    {
      ::glTranslated(0,0,-0.9);
      for(unsigned int ie=0;ie<cad.topo.aEdge.size();++ie){
        unsigned int iv0 = cad.topo.aEdge[ie].iv0;
        unsigned int iv1 = cad.topo.aEdge[ie].iv1;
        dfm2::CVec2d p = (cad.aVtx[iv0].pos+cad.aVtx[iv1].pos)*0.5;
        glyph.DrawStringAt(std::to_string(ie),0.001, p.x, p.y);
      }
      ::glTranslated(0,0,+0.9);
    }
    glfwSwapBuffers(viewer2.window);

    glfwPollEvents();
    viewer3.ExitIfClosed();
  }
}
