/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdlib.h>
#include <math.h>
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

#include "delfem2/v23m3q.h"
#include "delfem2/objfunc_v23.h"
#include "delfem2/objfunc_v23dtri.h"
#include "delfem2/v23m3q.h"
#include "delfem2/dtri_v2.h"
#include "delfem2/cad2d.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/rig_v3q.h"

#include "delfem2/../../external/tinygltf/tiny_gltf.h"
#include "delfem2/../../external/io_gltf.h"

// ----------------------------

#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_v23dtricad.h"
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl_rig_v23q.h"
#include "../glut_cam.h"

namespace dfm2 = delfem2;

// --------------------------------------------

std::vector<dfm2::CEPo2> aPo2D;
std::vector<dfm2::ETri> aETri;
std::vector<CVector2> aVec2;
std::vector<unsigned int> aLine;
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aXYZt;
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
std::vector<CInfoNearest> aInfoNearest;

std::vector<double> aXYZ_Contact;
std::vector<unsigned int> aTri_Contact;
std::vector<double> aNorm_Contact(aXYZ.size());
CBVH_MeshTri3D<dfm2::CBV3D_Sphere> bvh;
std::vector<double> aXYZ0_Contact;
std::vector<double> aRigWeight_Contact;
std::vector<unsigned int> aRigJoint_Contact;
std::vector<CRigBone> aBone;

const double dt = 0.01;
const double gravity[3] = {0.0, 0.0, 0.0};
const double contact_clearance = 0.0001;
const double rad_explore = 0.1;

CNav3D_GLUT nav;
bool is_animation = false;

// ------------------------------------

void StepTime()
{
  PBD_Pre3D(aXYZt,
            dt, gravity, aXYZ, aUVW, aBCFlag);
  PBD_TriStrain(aXYZt.data(),
                aXYZt.size()/3, aETri, aVec2);
  PBD_Bend(aXYZt.data(),
           aXYZt.size()/3, aETri, aVec2);
  PBD_Seam(aXYZt.data(),
           aXYZt.size()/3, aLine.data(), aLine.size()/2);
  Project_PointsIncludedInBVH_Outside_Cache(aXYZt.data(),aInfoNearest,
                                            aXYZt.size()/3,
                                            contact_clearance,bvh,
                                            aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                                            aTri_Contact.data(), aTri_Contact.size()/3,
                                            aNorm_Contact.data(), rad_explore);
  PBD_Post(aXYZ, aUVW,
           dt, aXYZt, aBCFlag);

}

// ---------------------------------------

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  ::glutPostRedisplay();
}

void myGlutDisplay(void)
{
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  //  ::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  nav.SetGL_Camera();
  
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

  
  ShowFPS();
  
  glutSwapBuffers();
}


void myGlutIdle(){
  if( is_animation ){
    StepTime();
  }
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  nav.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button, state, x,y);
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
    case 'q':
    case 'Q':
    case '\033':  /* '\033' ÇÕ ESC ÇÃ ASCII ÉRÅ[Éh */
      exit(0);
      break;
    case 'a':
    {
      is_animation = !is_animation;
      break;
    }
  }
}

void myGlutSpecial(int key, int x, int y){
  nav.glutSpecial(key,x,y);
}


class CRigidTrans_2DTo3D
{
public:
  CVector2 org2;
  CVector3 org3;
  CMatrix3 R;
};

int main(int argc,char* argv[])
{
  // Initailze GLUT
  ::glutInitWindowPosition(200,200);
  ::glutInitWindowSize(400, 300);
  ::glutInit(&argc, argv);
  ::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  ::glutCreateWindow("Cad View");
  
  // Set callback function
  ::glutMotionFunc(myGlutMotion);
  ::glutMouseFunc(myGlutMouse);
  ::glutDisplayFunc(myGlutDisplay);
  ::glutReshapeFunc(myGlutResize);
  ::glutKeyboardFunc(myGlutKeyboard);
  ::glutSpecialFunc(myGlutSpecial);
  ::glutIdleFunc(myGlutIdle);
  
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
    rt23.org2 = CVector2(2.5,0.5);
    rt23.org3 = CVector3(0.0, 0.0, 0.5);
    rt23.R.SetRotMatrix_Cartesian(0.0, 3.1415, 0.0);
    std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,1), cad);
    for(int iip=0;iip<aIP.size();++iip){
      const int ip = aIP[iip];
      CVector3 p0(aVec2[ip].x-rt23.org2.x, aVec2[ip].y-rt23.org2.y,0.0);
      CVector3 p1 = rt23.org3+rt23.R*p0;
      aXYZ[ip*3+0] = p1.x;
      aXYZ[ip*3+1] = p1.y;
      aXYZ[ip*3+2] = p1.z;
    }
    {
      CRigidTrans_2DTo3D rt23;
      rt23.org2 = CVector2(0.5,0.5);
      rt23.org3 = CVector3(0.0, 0.0, -0.5);
      rt23.R.SetIdentity();
      std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,0), cad);
      for(int iip=0;iip<aIP.size();++iip){
        const int ip = aIP[iip];
        CVector3 p0(aVec2[ip].x-rt23.org2.x, aVec2[ip].y-rt23.org2.y,0.0);
        CVector3 p1 = rt23.org3+rt23.R*p0;
        aXYZ[ip*3+0] = p1.x;
        aXYZ[ip*3+1] = p1.y;
        aXYZ[ip*3+2] = p1.z;
      }
    }
    aLine.clear();
    {
      std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(1, true, cad);
      std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(7, true, cad);
      const int npe = aIP0.size();
      assert( aIP1.size() == npe );
      for(int iip=0;iip<npe;++iip){
        int ip0 = aIP0[iip];
        int ip1 = aIP1[npe-iip-1];
        aLine.push_back(ip0);
        aLine.push_back(ip1);
      }
    }
    {
      std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(3, true, cad);
      std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(5, true, cad);
      const int npe = aIP0.size();
      assert( aIP1.size() == npe );
      for(int iip=0;iip<npe;++iip){
        int ip0 = aIP0[iip];
        int ip1 = aIP1[npe-iip-1];
        aLine.push_back(ip0);
        aLine.push_back(ip1);
      }
    }
  }
  aXYZt = aXYZ;
  
  { // make a unit sphere
//    MeshTri3D_Sphere(aXYZ_Contact, aTri_Contact, 0.3, 32, 32);
//    Rotate(aXYZ_Contact, 0.2, 0.3, 0.4);
    
    {
      tinygltf::Model model;
      tinygltf::TinyGLTF loader;
      std::string err;
      std::string warn;
      
      //std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.gltf";
      //bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, path_gltf);
      
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
      //      std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Monster.glb";
      std::string path_gltf = std::string(PATH_INPUT_DIR)+"/CesiumMan.glb";
      
      bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, path_gltf); // for binary glTF(.glb)
      if (!warn.empty()) { printf("Warn: %s\n", warn.c_str()); }
      if (!err.empty()) { printf("Err: %s\n", err.c_str()); }
      if (!ret) { printf("Failed to parse glTF\n"); return -1; }
      
      Print(model);
      ////
      GetMeshInfo(aXYZ0_Contact, aTri_Contact, aRigWeight_Contact, aRigJoint_Contact,
                  model, 0, 0);
      if( !model.skins.empty() ){
        aBone.resize( model.skins[0].joints.size() );
        std::vector<int> mapNode2Bone( model.nodes.size(), -1);
        for(int ij=0;ij<model.skins[0].joints.size();++ij){
          int inode = model.skins[0].joints[ij];
          mapNode2Bone[inode] = ij;
        }
        SetBone(aBone,
                model, model.skins[0].skeleton, -1, mapNode2Bone);
        GetBoneBinding(aBone,
                       model);
      }
      {
        Quat_Bryant(aBone[0].rot, -90.0/180.0*3.1415, 0, 0);
        aBone[0].trans[2] -= 1.0;
//        Quat_Bryant(aBone[0].rot, 0, 0, 0);
      }
      UpdateBoneRotTrans(aBone);
      aXYZ_Contact = aXYZ0_Contact;
      UpdateRigSkin(aXYZ_Contact.data(),
                    aXYZ0_Contact.data(), aXYZ0_Contact.size()/3,
                    aTri_Contact.data(), aTri_Contact.size()/3,
                    aBone, aRigWeight_Contact.data(), aRigJoint_Contact.data());
    }
    aNorm_Contact.resize(aXYZ_Contact.size());
    Normal_MeshTri3D(aNorm_Contact.data(),
                     aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                     aTri_Contact.data(), aTri_Contact.size()/3);
    bvh.Init(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
             aTri_Contact.data(), aTri_Contact.size()/3,
             0.01);
  }

  
  nav.camera.view_height = 1.0;
  nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  
  delfem2::opengl::setSomeLighting();
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
