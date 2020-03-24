/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <set>
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshtopo.h"
#include "delfem2/dtri.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/color.h"
//
#include "delfem2/v23m3q.h"
#include "delfem2/objfunc_v23.h"
#include "delfem2/objfunc_v23dtri.h"
#include "delfem2/dtri_v2.h"
#include "delfem2/cad2d_v2dtri.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/rig_v3q.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"

// ----------------------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_v23dtricad.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_rig_v23q.h"
#include "delfem2/opengl/glold_color.h"

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

void myGlutDisplay(const std::vector<dfm2::CDynTri>& aETri,
                   const std::vector<dfm2::CVec2d>& aVec2)
{
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  //  ::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  ::glPointSize(5);
  ::glLineWidth(1);
  {
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.8, 0.8, 0.8);
    /*
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    ::glEnable(GL_DEPTH_TEST);
     */
//    DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  }
  

//  DrawSphere_Edge(rad0);
}

class CRigidTrans_2DTo3D
{
public:
  dfm2::CVec3d Transform(const dfm2::CVec2d& pi) const {
    dfm2::CVec3d p0(pi.x()-org2.x(), pi.y()-org2.y(),0.0);
    dfm2::CVec3d p2 = p0;
    if( radinv_x < 1.0e-5 ) {
      double x0 = p0.x();
      p2.p[0] = x0;
      p2.p[2] = -0.5*radinv_x*x0*x0;
    }
    else{
      double x0 = p0.x();
      p2.p[0] = (1.0/radinv_x)*sin(radinv_x*x0);
      p2.p[2] = -(1.0/radinv_x)*(1-cos(radinv_x*x0));
    }
    dfm2::CVec3d p3 = org3+ dfm2::MatVec(R,p2);
    return p3;
  }
public:
  CRigidTrans_2DTo3D(){
    radinv_x = 0.0;
    R.SetIdentity();
  }
  dfm2::CVec2d org2;
  dfm2::CVec3d org3;
  dfm2::CMat3d R;
  double radinv_x;
};

int main(int argc,char* argv[])
{
  // -----------------------------
  // below: input data
  delfem2::CCad2D cad;
  double scale_adjust = 1.7;
  {
    std::string path_svg = std::string(PATH_INPUT_DIR)+"/tshirt.svg";
    dfm2::ReadSVG_Cad2D(cad, path_svg, 0.001*scale_adjust);
  }
  std::vector<unsigned int> aIESeam = {
    15, 6,
    13, 0,
    4, 9,
    11, 2,
    20, 17,
    22, 25,
    1, 18,
    12, 19,
    5, 24,
    8, 23
  };
  delfem2::CMesher_Cad2D mesher;
  mesher.edge_length = 0.015;
  std::vector<CRigidTrans_2DTo3D> aRT23;
  { // initial position
    aRT23.resize(4);
    { // back body
      CRigidTrans_2DTo3D& rt23 = aRT23[0];
      rt23.org2 = dfm2::CVec2d(0.189,-0.5)*scale_adjust;
      rt23.org3 = dfm2::CVec3d(0.0, 0.1, -0.2);
      rt23.R.SetRotMatrix_Cartesian(0.0, 3.1415, 0.0);
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[1];
      rt23.org2 = dfm2::CVec2d(0.506,-0.5)*scale_adjust;
      rt23.org3 = dfm2::CVec3d(0.0, 0.1, +0.2);
      rt23.R.SetIdentity();
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[2];
      rt23.org2 = dfm2::CVec2d(0.833,-0.45)*scale_adjust;
      rt23.org3 = dfm2::CVec3d(+0.3, 0.3, +0.0);
      rt23.R.SetRotMatrix_BryantAngle(-M_PI*0.5, +M_PI*0.5, 0.0);
      rt23.radinv_x = 13;
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[3];
      rt23.org2 = dfm2::CVec2d(1.148,-0.45)*scale_adjust;
      rt23.org3 = dfm2::CVec3d(-0.3, 0.3, +0.0);
      rt23.R.SetRotMatrix_BryantAngle(-M_PI*0.5, -M_PI*0.5, 0.0);
      rt23.radinv_x = 13;
    }

  }
  // above: input data
  // -----------------------------------
  // below: data preparation (derived)
  
  std::vector<dfm2::CDynTri> aETri;
  std::vector<dfm2::CVec2d> aVec2;
  { // make the seam edge equal number of division
    const double el = mesher.edge_length;
    for(int ie=0;ie<aIESeam.size()/2;++ie){
      const unsigned int ie0 = aIESeam[ie*2+0];
      const unsigned int ie1 = aIESeam[ie*2+1];
      const double l0 = cad.aEdge[ie0].LengthMesh();
      const double l1 = cad.aEdge[ie1].LengthMesh();
      const unsigned int ndiv = (int)(0.5*(l0+l1)/el+1);
      mesher.mapIdEd_NDiv[ie0] = ndiv;
      mesher.mapIdEd_NDiv[ie1] = ndiv;
    }
    delfem2::CMeshDynTri2D dmesh;
    mesher.Meshing(dmesh,
                   cad);
    dmesh.Check();
//    std::vector<dfm2::CDynPntSur> aPo2D;
//    aPo2D = dmesh.aEPo;
    aETri = dmesh.aETri;
    aVec2 = dmesh.aVec2;
  }
  std::vector<double> aXYZ; // deformed vertex positions
  std::vector<double> aXYZt;
  std::vector<double> aUVW; // deformed vertex velocity
  std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
  std::vector<dfm2::CInfoNearest<double>> aInfoNearest;
  {
    const int np = aVec2.size();
    aUVW.resize(np*3,0.0);
    aBCFlag.resize(np,0);
    aXYZ.resize(np*3);
  }
  for(int ifc=0;ifc<cad.aFace.size();++ifc){
    const CRigidTrans_2DTo3D& rt23 = aRT23[ifc];
    std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,ifc), cad);
    for(int ip : aIP){
      rt23.Transform(aVec2[ip]).CopyTo(aXYZ.data()+ip*3);
    }
  }
  std::vector<unsigned int> aLine;
  for(int ie=0;ie<aIESeam.size()/2;++ie){
    unsigned int ie0 = aIESeam[ie*2+0];
    unsigned int ie1 = aIESeam[ie*2+1];
    std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(ie0, true, cad);
    std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(ie1, true, cad);
    const int npe = aIP0.size();
    assert( aIP1.size() == npe );
    for(int iip=0;iip<npe;++iip){
      int ip0 = aIP0[iip];
      int ip1 = aIP1[npe-iip-1];
      aLine.push_back(ip0);
      aLine.push_back(ip1);
    }
  }
  aXYZt = aXYZ;
  
  // ----------
  
  std::vector<double> aXYZ_Contact;
  std::vector<unsigned int> aTri_Contact;
  std::vector<double> aNorm_Contact(aXYZ.size());
  dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere,double> bvh;
  std::vector<double> aXYZ0_Contact;
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
    aXYZ_Contact = aXYZ0_Contact;
    aNorm_Contact.resize(aXYZ_Contact.size());
    delfem2::Normal_MeshTri3D(aNorm_Contact.data(),
                              aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                              aTri_Contact.data(), aTri_Contact.size()/3);
    bvh.Init(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
             aTri_Contact.data(), aTri_Contact.size()/3,
             0.01);
  }
  
  // above: data preparation (derived)
  // ----------------------------------------------
  // below: opengl and UI

  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  // Enter main loop
  
  while (true)
  {
    StepTime(aXYZ, aXYZt, aUVW, aBCFlag, aInfoNearest,
             aETri,aVec2,aLine,
             aXYZ_Contact,aTri_Contact,aNorm_Contact,bvh);
    // ------------
    viewer.DrawBegin_oldGL();
//    myGlutDisplay(aETri,aVec2);
    for( auto& rt : aRT23 ){
      ::glPointSize(10);
      ::glColor3d(0,0,0);
      ::glBegin(GL_POINTS);
      dfm2::CVec3d v = rt.org3;
      ::glVertex3dv(v.p);
      ::glEnd();
    }
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
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
