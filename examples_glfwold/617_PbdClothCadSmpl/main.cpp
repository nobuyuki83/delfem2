/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <set>
#include "delfem2/garment.h"
#include "delfem2/objf_geo3.h"
#include "delfem2/geo3_v23m34q.h"

#include "delfem2/mshmisc.h"
#include "delfem2/mshtopo.h"
#include "delfem2/dtri.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/color.h"
//
#include "delfem2/objfdtri_objfdtri23.h"
#include "delfem2/cad2_dtri2.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/rig_geo3.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"
#include "inputs_garment.h"

// ----------------------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/caddtri_v3_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/rigv3_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

#ifndef M_PI
#  define  M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// --------------------------------------------
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

// ---------------------------------------

/*
namespace cereal{

template<class Archive>
void serialize(Archive & archive, dfm2::CRigidTrans_2DTo3D &rt23)
{
  archive( cereal::make_nvp("radinv_x", rt23.radinv_x) );
  archive( cereal::make_nvp("org2",     rt23.org2.p) );
  archive( cereal::make_nvp("org3",     rt23.org3.p) );
  archive( cereal::make_nvp("R",        rt23.R.mat) );
}

}
 */





void Draw
(const std::vector<dfm2::CDynTri>& aETri,
 const std::vector<double>& aXYZ,
 std::vector<double>& aXYZ_Contact,
 std::vector<unsigned int>& aTri_Contact)
{
  //    myGlutDisplay(aETri,aVec2);
  /*
  ::glDisable(GL_LIGHTING);
  for( auto& rt : aRT23 ){
    ::glPointSize(10);
    ::glColor3d(0,0,0);
    ::glBegin(GL_POINTS);
    dfm2::CVec3d v = rt.org3;
    ::glVertex3dv(v.p);
    ::glEnd();
  }
   */
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


int main(int argc,char* argv[])
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
  
  
  // above: data preparation (derived)
  // ----------------------------------------------
  // below: opengl and UI

  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  while (true)
  {
    StepTime(aXYZ, aXYZt, aUVW, aBCFlag, aInfoNearest_Contact,
             aETri,aVec2,aLine,
             aXYZ_Contact,aTri_Contact,aNorm_Contact,bvh_Contact);
    // ------------
    viewer.DrawBegin_oldGL();
    Draw(aETri,aXYZ,
         aXYZ_Contact,aTri_Contact);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
