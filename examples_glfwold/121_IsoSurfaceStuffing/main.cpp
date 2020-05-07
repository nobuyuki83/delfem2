/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <ctime>
#include "delfem2/primitive.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/iss.h"
#include "delfem2/bv.h"
#include "delfem2/srch_v3bvhmshtopo.h"

// ---------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ------------------------------------------


std::vector<double> aXYZ;
std::vector<unsigned int> aTet;
std::vector<unsigned int> aTetSurface;
std::vector<dfm2::CColor> aTetColor;

std::vector<unsigned int> aTet1;
std::vector<dfm2::CColor> aTetColor1;

double vis_cut_org[3] = {-0.0, 0.0, 0.0};
double vis_cut_nrm[3] = {0.0,-0.9, +0.2};
bool is_animation = false;
double cur_time = 0.5;
int imode_draw = 0;


void SetProblem(int iprob)
{
  std::vector<int> aIsOnSurfXYZ;
  if( iprob == 0 ){
    class CInSphere : public dfm2::CInput_IsosurfaceStuffing
    {
    public:
      explicit CInSphere(double rad){
        sp.radius_ = rad;
      }
      double SignedDistance(double px, double py, double pz) const override {
        double n[3];
        return sp.Projection(n,
                             px,py,pz);
      }
      void Level(
          int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
          double px, double py, double pz) const override
      {
        sdf = this->SignedDistance(px,py,pz);
        const double radius0 = sp.radius_;
        const double radius2 = sqrt(px*px+py*py+pz*pz);
        if( radius2 > radius0*0.5 ){ ilevel_vol = 0; }
        else{ ilevel_vol = 1; }
        ilevel_srf = 1;
        nlayer = 1;
        ilevel_vol = -1;
      }
    public:
      dfm2::CSphere<double> sp;
    };
    double rad = 1.5;
    CInSphere sphere(rad);
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                       sphere, 0.2, rad*4.0, cent);
  }
  else if( iprob ==  1 ){
    class CInBox : public dfm2::CInput_IsosurfaceStuffing
    {
    public:
      CInBox(double hwx, double hwy, double hwz){
        bx.hwx = hwx;
        bx.hwy = hwy;
        bx.hwz = hwz;
      }
      double SignedDistance(double px, double py, double pz) const override {
        double n[3];
        return bx.Projection(n,
                             px,py,pz);
      }
      void Level(
          int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
          double px, double py, double pz) const override
      {
        sdf = this->SignedDistance(px,py,pz);
        ilevel_vol = -1;
        ilevel_srf = 1;
        nlayer = 1;
      }
    public:
      dfm2::CBox<double> bx;
    };
    const double hwx = 0.91;
    const double hwy = 0.61;
    const double hwz = 0.41;
    double cent[3] = {0,0,0};
    CInBox box(hwx,hwy,hwz);
    dfm2::IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                       box, 0.1, 2.00, cent);
  }
  else if( iprob == 2 ){
    class CCavSphere : public dfm2::CInput_IsosurfaceStuffing
    {
    public:
      CCavSphere(){
        const double hwx = 0.91;
        const double hwy = 0.61;
        const double hwz = 0.61;
        box.hwx = hwx;
        box.hwy = hwy;
        box.hwz = hwz;
        ////
        const double rad = 0.2;
        sphere.radius_ = rad;
      }
      double SignedDistance(double x, double y, double z ) const override {
        double n[3];
        double dist0 = -sphere.Projection(n,
                                          x, y, z);
        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;
        double dist1 = box.Projection(n,
                                      x-cx, y-cy, z-cz);
        return (dist0<dist1) ? dist0 : dist1;
      }
      void Level(
          int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
          double px, double py, double pz) const override
      {
        sdf = this->SignedDistance(px,py,pz);
        double dist0 = sqrt(px*px+py*py+pz*pz);
        ilevel_vol = -1;
        if( dist0 < 0.5 && px > 0 ){ ilevel_srf = +3; }
        else{                        ilevel_srf = -1; }
        nlayer = 1;
      }
    public:
      dfm2::CBox<double> box;
      dfm2::CSphere<double> sphere;
    } cav_sphere;
    double cent[3] = {0,0,0};
    dfm2::IsoSurfaceStuffing(aXYZ, aTet, aIsOnSurfXYZ,
                       cav_sphere, 0.2, 3.00, cent);
  }
  if( iprob == 3 ){
    class CMesh : public dfm2::CInput_IsosurfaceStuffing
    {
    public:
      double SignedDistance(double x, double y, double z) const override {
        dfm2::CVec3d n0;
        return obj.SignedDistanceFunction(n0,
                                          dfm2::CVec3d(x,y,z), aXYZ_Tri, aTri, aNorm);
      }
      void Level(
          int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
          double px, double py, double pz) const override
      {
        sdf = this->SignedDistance(px,py,pz);
        ilevel_vol = 0;
        ilevel_srf = 2;
        nlayer = 2;
      }
    public:
      dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> obj;
      std::vector<double> aNorm;
      std::vector<unsigned int> aTri;
      std::vector<double> aXYZ_Tri;
    } mesh;
    {
      delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
               mesh.aXYZ_Tri, mesh.aTri);
      delfem2::Normalize_Points3(mesh.aXYZ_Tri,2.3);
      mesh.obj.Init(mesh.aXYZ_Tri.data(), mesh.aXYZ_Tri.size()/3,
                    mesh.aTri.data(), mesh.aTri.size()/3,
                    0.0);
      mesh.aNorm.resize(mesh.aXYZ_Tri.size());
      delfem2::Normal_MeshTri3D(mesh.aNorm.data(),
                                mesh.aXYZ_Tri.data(), mesh.aXYZ_Tri.size()/3,
                                mesh.aTri.data(), mesh.aTri.size()/3);
    }
    double cent[3] = {0,0,0};
    dfm2::IsoSurfaceStuffing(aXYZ, aTet, aIsOnSurfXYZ,
                             mesh, 0.18, 3.0, cent);
  }
  
  aTetColor.resize(aTet.size()/4);
  for(size_t it=0;it<aTet.size()/4;++it){
    aTetColor[it].setRandomVividColor();
  }
  
  std::vector<int> aTetSurRel;
  ElSuEl_MeshElem(aTetSurRel,
                              aTet.data(), aTet.size()/4,
                              delfem2::MESHELEM_TET,
                              aXYZ.size()/3);
  aTetSurface.clear();
  for(size_t it=0;it<aTet.size()/4;++it){
    for(int ift=0;ift<4;++ift){
      if( aTetSurRel[it*8+ift*2+0] != -1 ){ continue; }
      aTetSurface.push_back(it);
      aTetSurface.push_back(ift);
    }
  }
  
  aTet1.clear();
  aTetColor1.clear();
}

void myGlutDisplay()
{
  ::glEnable(GL_LIGHTING);
  if( imode_draw == 0 ){
    dfm2::opengl::DrawMeshTet3D_Cut(aXYZ,aTet,aTetColor,
                              vis_cut_org, vis_cut_nrm);
  }
  else if( imode_draw == 1 ){
    dfm2::opengl::DrawMeshTet3DSurface_Edge(aXYZ, aTet, aTetSurface);
  }
  else if( imode_draw == 2 ){
    dfm2::opengl::DrawMeshTet3D_Cut(aXYZ,aTet1,aTetColor1,
                                    vis_cut_org, vis_cut_nrm);
  }
}

// ------------------------------------

int main(int argc,char* argv[])
{	
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
	
  viewer.nav.camera.view_height = 2;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  ////
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_LIGHT0);            
  float light0pos[4] = {0.5,0.5,+20,0};
  ::glLightfv(GL_LIGHT0, GL_POSITION, light0pos);      
  float white[3] = {1.0,1.0,1.0};
  ::glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
  
  ::glEnable(GL_LIGHT1);    
  float light1pos[4] = {0.5,0.5,-20,0}; 
  ::glLightfv(GL_LIGHT1, GL_POSITION, light1pos);          
  float white2[3] = {1.0,1.0,0.0};
  ::glLightfv(GL_LIGHT1, GL_DIFFUSE, white2);        
  ///
  
  while(!glfwWindowShouldClose(viewer.window)){
    {
      static int iframe = 0;
      int iprob = iframe/100;
      if( iframe % 100 == 0 ){
        SetProblem(iprob);
      }
      iframe = (iframe+1)%300;
      //
      cur_time += 0.02;
      if( cur_time > 1 ){ cur_time = 0.0; }
      vis_cut_org[1] = -1*(1-cur_time) + 1*cur_time;
    }
    // ------------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

