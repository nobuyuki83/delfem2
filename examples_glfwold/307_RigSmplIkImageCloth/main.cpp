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
#include <random>
#include <ctime>
#include <GLFW/glfw3.h>
#include "delfem2/rigopt.h"
#include "delfem2/garment.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/mshmisc.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/objf_v23m34q.h"
#include "delfem2/objf_dtriv23m34q.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/caddtri_v3_glold.h"
#include "delfem2/cnpy/smpl_cnpy.h"
#include "inputs_garment.h"

#include "delfem2/rig_v3q.h"
#include "delfem2/opengl/rigv3_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/tex_gl.h"

#include "inputs_imgboneloc.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace dfm2 = delfem2;

// -------------------

const double dt = 0.01;
const double gravity[3] = {0.0, -0.1, 0.0};
const double contact_clearance = 0.0001;
const double rad_explore = 0.1;

// ------------------------------------

void StepTime
(std::vector<double>& aXYZ, // deformed vertex positions
 std::vector<double>& aXYZt,
 std::vector<double>& aUVW, // deformed vertex velocity
 std::vector<int>& aBCFlag,  // boundary condition flag (0:free 1:fixed)
 std::vector<dfm2::CInfoNearest<double>>& aInfoNearest,
 const std::vector<dfm2::CDynTri>& aETri,
 const std::vector<dfm2::CVec2d>& aVec2,
 const std::vector<unsigned int>& aLine,
 //
 std::vector<double>& aXYZ_Contact,
 std::vector<unsigned int>& aTri_Contact,
 std::vector<double>& aNorm_Contact,
 dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere,double>& bvh)
{
  dfm2::PBD_Pre3D(aXYZt,
                  dt, gravity, aXYZ, aUVW, aBCFlag);
  dfm2::PBD_TriStrain(aXYZt.data(),
                      aXYZt.size()/3, aETri, aVec2);
  dfm2::PBD_Bend(aXYZt.data(),
                 aXYZt.size()/3, aETri, aVec2, 0.5);
  dfm2::PBD_Seam(aXYZt.data(),
                 aXYZt.size()/3, aLine.data(), aLine.size()/2);
  dfm2::Project_PointsIncludedInBVH_Outside_Cache(aXYZt.data(),aInfoNearest,
                                                  aXYZt.size()/3,
                                                  contact_clearance,bvh,
                                                  aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                                                  aTri_Contact.data(), aTri_Contact.size()/3,
                                                  aNorm_Contact.data(), rad_explore);
  dfm2::PBD_Post(aXYZ, aUVW,
                 dt, aXYZt, aBCFlag);
}


void Draw
(const std::vector<double>& aXYZ1,
 const std::vector<unsigned int>& aTri,
 const std::vector<dfm2::CRigBone>& aBone,
 const std::vector<dfm2::CTarget>& aTarget)
{
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_DEPTH_TEST);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
  //    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
  ::glDisable(GL_DEPTH_TEST);
  ::glDisable(GL_LIGHTING);
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for(int it=0;it<aTarget.size();++it){
    const unsigned int ib = aTarget[it].ib;
    ::glColor3d(1,0,0);
    dfm2::opengl::myGlVertex(aBone[ib].Pos());
  }
  ::glEnd();
  // ------
  ::glEnable(GL_DEPTH_TEST);
  ::glBegin(GL_LINES);
  ::glColor3d(1,0,0);
  for(int it=0;it<aTarget.size();++it){
    dfm2::CVec3d p = aTarget[it].pos;
    dfm2::opengl::myGlVertex(p+10*dfm2::CVec3d(0,0,1));
    dfm2::opengl::myGlVertex(p-10*dfm2::CVec3d(0,0,1));
  }
  ::glEnd();
  /*
   ::glDisable(GL_DEPTH_TEST);
   delfem2::opengl::DrawBone(aBone,
   -1, -1,
   0.01, 1.0);
   */
  ::glEnable(GL_DEPTH_TEST);
  //    dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat1,0.02);
}



// --------------------

int main()
{
  // -----------------------------
  // below: input data
  std::vector<dfm2::CDynTri> aETri;
  std::vector<dfm2::CVec2d> aVec2;
  std::vector<double> aXYZ; // deformed vertex positions
  std::vector<unsigned int> aLine;
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
    dfm2::MeshingPattern(aETri,aVec2,aXYZ,aLine,
                         aRT23,cad,aIESeam,mesher_edge_length);
  }
  std::vector<double> aXYZt = aXYZ;
  std::vector<double> aUVW(aXYZ.size(), 0.0);
  std::vector<int> aBCFlag(aXYZ.size()/3, 0.0);
  
  // ----------
  
  std::vector<double> aXYZ0_Contact;
  std::vector<unsigned int> aTri_Contact;
  {
    std::vector<dfm2::CRigBone> aBone;
    { // makineg aBone
      std::vector<int> aIndBoneParent;
      std::vector<double> aJntRgrs0;
      std::vector<double> aRigWeight_Contact;
      dfm2::cnpy::LoadSmpl(aXYZ0_Contact,
                           aRigWeight_Contact,
                           aTri_Contact,
                           aIndBoneParent,
                           aJntRgrs0,
                           std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
      dfm2::Smpl2Rig(aBone,
                     aIndBoneParent, aXYZ0_Contact, aJntRgrs0);
      dfm2::UpdateBoneRotTrans(aBone);
    }
  }
  std::vector<double> aXYZ_Contact = aXYZ0_Contact;
  std::vector<double> aNorm_Contact(aXYZ_Contact.size());
  delfem2::Normal_MeshTri3D(aNorm_Contact.data(),
                            aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                            aTri_Contact.data(), aTri_Contact.size()/3);
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere,double> bvh_Contact;
  bvh_Contact.Init(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                   aTri_Contact.data(), aTri_Contact.size()/3,
                   0.01);
  std::vector<dfm2::CInfoNearest<double>> aInfoNearest_Contact;
  
  std::vector<dfm2::CTarget> aTarget;
  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    std::string name_img_in_test_inputs;
    std::vector< std::pair<double,dfm2::CVec2d> > aBoneLoc;
    double scale = 1.0;
    BoneLocs_SmplUglysweater(name_img_in_test_inputs,
                             scale,
                             aBoneLoc);
    //--
    int width, height;
    {
      int channels;
      unsigned char *img = stbi_load((std::string(PATH_INPUT_DIR)+"/"+name_img_in_test_inputs).c_str(),
                                     &width, &height, &channels, 0);
      tex.Initialize(width, height, img, "rgb");
      delete[] img;
      tex.max_x = -scale*width*0.5;
      tex.min_x = +scale*width*0.5;
      tex.max_y = -scale*height*0.5;
      tex.min_y = +scale*height*0.5;
      tex.z = -0.5;
    }
    aTarget.clear();
    for(int it=0;it<aBoneLoc.size();++it){
      dfm2::CTarget t;
      t.ib = aBoneLoc[it].first;
      int iw = aBoneLoc[it].second.x();
      int ih = aBoneLoc[it].second.y();
      t.pos.p[0] = (double)iw/width-0.5;
      t.pos.p[1] = 0.5*height/width - (double)ih/height;
      aTarget.push_back(t);
    }
  }
  
  std::vector<double> aXYZ0_Body;
  std::vector<double> aW;
  std::vector<dfm2::CRigBone> aBone;
  std::vector<unsigned int> aTri;
  {
    std::vector<int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl(aXYZ0_Body,
                         aW,
                         aTri,
                         aIndBoneParent,
                         aJntRgrs,
                         std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    dfm2::Smpl2Rig(aBone,
                   aIndBoneParent, aXYZ0_Body, aJntRgrs);
    
  }
  
  std::vector<double> aXYZ1 = aXYZ0_Body;
  { // initalize pose
    dfm2::UpdateBoneRotTrans(aBone);
    dfm2::Skinning_LBS(aXYZ1,
                       aXYZ0_Body, aBone, aW);
  }
  std::vector< std::pair<dfm2::CVec3d,dfm2::CVec3d> > aTargetOriginPos;
  for(int it=0;it<aTarget.size();++it){
    unsigned int ib = aTarget[it].ib;
    aTargetOriginPos.push_back( std::make_pair(aTarget[it].pos,
                                               aBone[ib].Pos()) );
  }
  
  // -----------
//  Check(aBone, aTarget);
     
  // -----------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.camera_rot_mode = dfm2::CAMERA_ROT_YTOP;
  viewer.nav.camera.view_height = 1.0;
  dfm2::opengl::setSomeLighting();
  tex.InitGL();

  int iframe = 0;
  while (true)
  {
    {
      double r = iframe*0.01;
      if( r > 1 ){ r = 1; }
      for(int it=0;it<aTarget.size();++it){
        aTarget[it].pos = r*aTargetOriginPos[it].first + (1-r)*aTargetOriginPos[it].second;
      }
      Solve_MinRigging(aBone, aTarget);
      Skinning_LBS(aXYZ1,
                   aXYZ0_Body, aBone, aW);
    }
    if( iframe > 200 ){
      for(int ib=0;ib<aBone.size();++ib){
        dfm2::Quat_Identity(aBone[ib].quatRelativeRot);
      }
      dfm2::UpdateBoneRotTrans(aBone);
      for(int it=0;it<aTarget.size();++it){
        aTarget[it].pos = aTargetOriginPos[it].second;
      }
      Solve_MinRigging(aBone, aTarget);
      Skinning_LBS(aXYZ1,
                   aXYZ0_Body, aBone, aW);
      iframe = 0;
    }
    
    // -------------------
    viewer.DrawBegin_oldGL();
    Draw(aXYZ1,aTri,aBone,aTarget);
    tex.Draw_oldGL();
    {
      ::glEnable(GL_LIGHTING);
      dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Gray(0.8) );
      ::glColor3d(1,0,0);
      //    delfem2::opengl::DrawMeshTri3D_Edge(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
      //                                        aTri_Contact.data(), aTri_Contact.size()/3);
      ::glEnable(GL_LIGHTING);
      delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ_Contact.data(),
                                              aTri_Contact.data(), aTri_Contact.size()/3);
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0,0,0);
      //    delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ, aETri);
      delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ, aETri);
      ::glEnable(GL_LIGHTING);
      dfm2::opengl::myGlColorDiffuse( dfm2::CColor::Red() );
      delfem2::opengl::DrawMeshDynTri_FaceNorm(aETri, aXYZ.data());
    }
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
    iframe++;
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
