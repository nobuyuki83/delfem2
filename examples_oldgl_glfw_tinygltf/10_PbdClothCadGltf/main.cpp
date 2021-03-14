/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/tinygltf/io_gltf.h"
#include "tinygltf/tiny_gltf.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/opengl/glfw/viewer3.h"
#include "delfem2/opengl/old/caddtri_v3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/pbd_geo3.h"
#include "delfem2/pbd_geo3dtri23.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/cad2_dtri2.h"
#include "delfem2/dtri.h"
#include "delfem2/mshmisc.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/srchbv3sphere.h"
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <vector>
#include <set>

namespace dfm2 = delfem2;

// --------------------------------------------

std::vector<dfm2::CDynPntSur> aPo2D;
std::vector<dfm2::CDynTri> aETri;
std::vector<dfm2::CVec2d> aVec2;
std::vector<unsigned int> aLine;
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aXYZt;
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
std::vector<dfm2::CInfoNearest<double>> aInfoNearest;

std::vector<double> aXYZ_Contact;
std::vector<unsigned int> aTri_Contact;
std::vector<double> aNorm_Contact(aXYZ.size());
dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere,double> bvh;
std::vector<double> aXYZ0_Contact;
std::vector<double> aRigSparseW_Contact;
std::vector<unsigned int> aRigSparseI_Contact;
std::vector<dfm2::CRigBone> aBone;

const double dt = 0.01;
const double gravity[3] = {0.0, 0.0, 0.0};
const double contact_clearance = 0.0001;
const double rad_explore = 0.1;

bool is_animation = false;

// ------------------------------------

void StepTime()
{
  dfm2::PBD_Pre3D(
      aXYZt,
      dt, gravity, aXYZ, aUVW, aBCFlag);
  dfm2::PBD_TriStrain(
      aXYZt.data(),
      aXYZt.size()/3, aETri, aVec2);
  dfm2::PBD_Bend(
      aXYZt.data(),
      aXYZt.size()/3, aETri, aVec2, 1.0);
  dfm2::PBD_Seam(
      aXYZt.data(),
      aXYZt.size()/3, aLine.data(), aLine.size()/2);
  dfm2::Project_PointsIncludedInBVH_Outside_Cache(
      aXYZt.data(),aInfoNearest,
      aXYZt.size()/3,
      contact_clearance,bvh,
      aXYZ_Contact.data(), aXYZ_Contact.size()/3,
      aTri_Contact.data(), aTri_Contact.size()/3,
      aNorm_Contact.data(), rad_explore);
  dfm2::PBD_Post(
      aXYZ, aUVW,
      dt, aXYZt, aBCFlag);
}

// ---------------------------------------

void myGlutDisplay()
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
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ, aETri);

  ::glColor3d(1,0,0);
  delfem2::opengl::DrawMeshTri3D_Edge(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                                      aTri_Contact.data(), aTri_Contact.size()/3);
//  DrawSphere_Edge(rad0);
}

class CRigidTrans_2DTo3D
{
public:
  dfm2::CVec2d org2;
  dfm2::CVec3d org3;
  dfm2::CMat3d R;
};

int main(int argc,char* argv[])
{
  delfem2::CCad2D cad;
  {
    const double xys0[8] = { -0.0,-0.0,  +1.0,-0.0, +1.0,+1.0, -0.0,+1.0, };
    const double xys1[8] = { +2.0,-0.0,  +3.0,-0.0, +3.0,+1.0, +2.0,+1.0, };
    cad.AddPolygon(std::vector<double>(xys0,xys0+8));
    cad.AddPolygon(std::vector<double>(xys1,xys1+8));
  }
  delfem2::CMesher_Cad2D mesher;
  mesher.edge_length = 0.04;
  delfem2::CMeshDynTri2D dmesh;
  mesher.Meshing(dmesh,
                 cad);
  dmesh.Check();
  aPo2D = dmesh.aEPo;
  aETri = dmesh.aETri;
  aVec2 = dmesh.aVec2;

  // -----------------------------
  const int np = aPo2D.size();
  aUVW.resize(np*3,0.0);
  aBCFlag.resize(np,0);
  aXYZ.resize(np*3);
  {
    CRigidTrans_2DTo3D rt23;
    rt23.org2 = dfm2::CVec2d(2.5,0.5);
    rt23.org3 = dfm2::CVec3d(0.0, 0.0, 0.5);
    rt23.R.SetRotMatrix_Cartesian(0.0, 3.1415, 0.0);
    std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,1), cad);
    for(int ip : aIP){
      dfm2::CVec3d p0(aVec2[ip].x()-rt23.org2.x(), aVec2[ip].y()-rt23.org2.y(),0.0);
      dfm2::CVec3d p1 = rt23.org3+ dfm2::MatVec(rt23.R,p0);
      aXYZ[ip*3+0] = p1.x();
      aXYZ[ip*3+1] = p1.y();
      aXYZ[ip*3+2] = p1.z();
    }
    {
      CRigidTrans_2DTo3D rt23;
      rt23.org2 = dfm2::CVec2d(0.5,0.5);
      rt23.org3 = dfm2::CVec3d(0.0, 0.0, -0.5);
      rt23.R.SetIdentity();
      std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,0), cad);
      for(int ip : aIP){
        dfm2::CVec3d p0(aVec2[ip].x()-rt23.org2.x(), aVec2[ip].y()-rt23.org2.y(),0.0);
        dfm2::CVec3d p1 = rt23.org3+dfm2::MatVec(rt23.R,p0);
        aXYZ[ip*3+0] = p1.x();
        aXYZ[ip*3+1] = p1.y();
        aXYZ[ip*3+2] = p1.z();
      }
    }
    aLine.clear();
    {
      std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(1, true, cad);
      std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(7, true, cad);
      const unsigned int npe = aIP0.size();
      assert( aIP1.size() == npe );
      for(unsigned int iip=0;iip<npe;++iip){
        int ip0 = aIP0[iip];
        int ip1 = aIP1[npe-iip-1];
        aLine.push_back(ip0);
        aLine.push_back(ip1);
      }
    }
    {
      std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(3, true, cad);
      std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(5, true, cad);
      const std::size_t npe = aIP0.size();
      assert( aIP1.size() == npe );
      for(unsigned int iip=0;iip<npe;++iip){
        unsigned int ip0 = aIP0[iip];
        unsigned int ip1 = aIP1[npe-iip-1];
        aLine.push_back(ip0);
        aLine.push_back(ip1);
      }
    }
  }
  aXYZt = aXYZ;
  
  { // make a unit sphere
    {
      tinygltf::Model model;
      tinygltf::TinyGLTF loader;
      std::string err;
      std::string warn;
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Monster.glb";
      std::string path_gltf = std::string(PATH_INPUT_DIR)+"/CesiumMan.glb";
      bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, path_gltf); // for binary glTF(.glb)
      if (!warn.empty()) { printf("Warn: %s\n", warn.c_str()); }
      if (!err.empty()) { printf("Err: %s\n", err.c_str()); }
      if (!ret) { printf("Failed to parse glTF\n"); return -1; }
      dfm2::Print(model);
      dfm2::GetMeshInfo(
          aXYZ0_Contact, aTri_Contact,
          aRigSparseW_Contact, aRigSparseI_Contact,
          model, 0, 0);
      dfm2::GetBone(aBone, model, 0);
      {
        dfm2::Quat_Bryant(aBone[0].quatRelativeRot, -90.0/180.0*3.1415, 0.0, 0.0);
        aBone[0].transRelative[2] -= 1.0;
//        Quat_Bryant(aBone[0].rot, 0, 0, 0);
      }
      UpdateBoneRotTrans(aBone);
      aXYZ_Contact = aXYZ0_Contact;
      dfm2::Skinning_LBS_LocalWeight(
          aXYZ_Contact.data(),
          aXYZ0_Contact.data(), aXYZ0_Contact.size()/3,
          aBone, aRigSparseW_Contact.data(), aRigSparseI_Contact.data());
    }
    aNorm_Contact.resize(aXYZ_Contact.size());
    delfem2::Normal_MeshTri3D(
        aNorm_Contact.data(),
        aXYZ_Contact.data(), aXYZ_Contact.size()/3,
        aTri_Contact.data(), aTri_Contact.size()/3);
    bvh.Init(
        aXYZ_Contact.data(), aXYZ_Contact.size()/3,
        aTri_Contact.data(), aTri_Contact.size()/3,
        0.01);
  }

  delfem2::opengl::CViewer3 viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  // Enter main loop
  
  while (true)
  {
    StepTime();
    // ------------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
