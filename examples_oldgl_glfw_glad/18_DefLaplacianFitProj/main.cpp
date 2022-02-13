/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should be before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <glad/glad.h> // this should come first
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/r2tglo.h"
#include "delfem2/opengl/r2t.h"
#include "delfem2/deflap.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_primitive.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

class CCapsuleRigged
{
public:
  // -----
  void Init(){
    dfm2::MeshTri3_Capsule(aXYZ0,aElm,
                           0.2, 1.0,
                           16, 5, 8);
    // -------
    aBone.clear();
    {
      dfm2::CRigBone b;
      b.ibone_parent = -1;
      b.invBindMat[7] = +0.5;
      b.transRelative[1] = -0.5;
      aBone.push_back(b);
    }
    {
      dfm2::CRigBone b;
      b.ibone_parent = 0;
      b.transRelative[1] = +0.5;
      aBone.push_back(b);
    }
    // -------
    {
      unsigned int np = aXYZ0.size()/3;
      unsigned int nb = aBone.size();
      aW.resize(np*nb);
      for(unsigned int ip=0;ip<aXYZ0.size()/3;++ip) {
        const double* p0 = aXYZ0.data()+ip*3;
        double w_tot = 0;
        for(unsigned int ib=0;ib<nb;++ib){
          double pb[3] = {
            -aBone[ib].invBindMat[3],
            -aBone[ib].invBindMat[7],
            -aBone[ib].invBindMat[11]};
          double len = dfm2::Distance3(p0,pb);
          double wb = 1.0/(len+1.0e-10);
          aW[ip*nb+ib] = wb;
          w_tot += wb;
        }
        for(unsigned int ib=0;ib<nb;++ib) {
          aW[ip*nb+ib] /= w_tot;
        }
      }
    }
    aXYZ1 = aXYZ0;
  }
  void Def(double angle)
  {
    dfm2::Quat_Bryant(aBone[1].quatRelativeRot,  0.,0.,angle);
    dfm2::UpdateBoneRotTrans(aBone);
    dfm2::Skinning_LBS(aXYZ1,
        aXYZ0, aBone, aW);
  }
public:
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aElm;
  std::vector<double> aXYZ1;
  std::vector<dfm2::CRigBone> aBone;
  std::vector<double> aW;
};

void Project(
    std::vector< std::pair<unsigned int, dfm2::CVec3d> >& aIdpNrm,
    std::vector<double>& aXYZ1,
    const dfm2::opengl::CRender2Tex_DrawOldGL_BOX& sampler)
{
  aIdpNrm.clear();
  for(unsigned int ip=0;ip<aXYZ1.size()/3;++ip){
    const dfm2::CVec3d& ps = dfm2::CVec3d(aXYZ1.data()+ip*3);
    dfm2::CVec3d pmin, nmin;
    double len_min = -1;
    for(const auto & smplr : sampler.aSampler){
      dfm2::CVec3d p0, n0;
      bool res = GetProjectedPoint(p0, n0, ps, smplr);
      if( !res ){ continue; }
      dfm2::CVec3d z0(0,0,-1);
      dfm2::CMat4d mMVP = dfm2::CMat4d(smplr.mat_projection) * dfm2::CMat4d(smplr.mat_modelview);
      dfm2::CMat4d mMVPinv = mMVP.Inverse();
      dfm2::CVec3d z1 = mMVPinv.MultVec3_Homography(z0.p);
      double ct = n0.dot(z1.normalized());
      if( ct <= 0.0 ){ continue; }
      if( (p0-ps).norm() > 0.1 ){ continue; }
      const double len = ((p0-ps).norm()+1.0e-5)/ct;
      if( len_min < 0 || len < len_min ){
        len_min = len;
        pmin = p0;
        nmin = n0;
      }
    }
    if( len_min < 0 ){ continue; }
    {
      aIdpNrm.emplace_back( ip, nmin );
      aXYZ1[ip*3+0] = pmin.x;
      aXYZ1[ip*3+1] = pmin.y;
      aXYZ1[ip*3+2] = pmin.z;
    }
  }
}


void Draw(
    const CCapsuleRigged& trg,
    const std::vector<double>& aXYZ1,
    const std::vector<unsigned int>& aTri)
{
  dfm2::opengl::DrawBackground( dfm2::CColor(0.2,0.7,0.7) );
  ::glEnable(GL_LIGHTING);
  ::glColor3d(1,1,1);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1,aTri);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glLineWidth(1);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ1,aTri);
  {
    glPointSize(3);
    ::glLineWidth(1);
    dfm2::opengl::DrawMeshTri3D_Edge(trg.aXYZ1, trg.aElm);
  }
}

void LaplacianLinear(
    std::vector<double>& aXYZ1,
    CCapsuleRigged& trg,
    dfm2::opengl::CRender2Tex_DrawOldGL_BOX& sampler,
    dfm2::glfw::CViewer3& viewer,
    const std::vector<double>& aXYZ0,
    const std::vector<unsigned int>& aTri)
{
  glfwSetWindowTitle(viewer.window,"Laplacian Linear");
  dfm2::CDef_LaplacianLinear def;
  def.Init(aXYZ0, aTri, true);
  def.aBCFlag.assign(aXYZ0.size(),0);
  def.max_itr = 1000;
  def.weight_nrm = 1;
  for(unsigned int iframe=0;iframe<30;iframe++){
    // deform target mesh
    trg.Def(0.5*sin(0.1*iframe));

    // sample the space
    for(auto& smplr: sampler.aSampler){
      smplr.InitGL(); // move the sampled image to a texture
      smplr.Start();
      dfm2::opengl::SetView(smplr);
      ::glDisable(GL_POLYGON_OFFSET_FILL ); // the depth will be jazzy without this
      ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(trg.aXYZ1,trg.aElm);
      smplr.End();
    }

    // deformation
    Project(def.aIdpNrm, aXYZ1, sampler);
    def.SetValueToPreconditioner();
    def.Deform(aXYZ1, aXYZ0);
    std::cout << iframe << " nitr:" << def.aConvHist.size() << " nconst:" << def.aIdpNrm.size() << " np:" << aXYZ0.size()/3 << std::endl;

    // drawing functions
    viewer.DrawBegin_oldGL();
    Draw(trg,aXYZ1,aTri);
    sampler.Draw();
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

void LaplacianDegenerate(
    std::vector<double>& aXYZ1,
    CCapsuleRigged& trg,
    dfm2::opengl::CRender2Tex_DrawOldGL_BOX& sampler,
    dfm2::glfw::CViewer3& viewer,
    const std::vector<double>& aXYZ0,
    const std::vector<unsigned int>& aTri)
{ // test degenerate deformer
  glfwSetWindowTitle(viewer.window,"Laplacian Degenerate");
  dfm2::CDef_LaplacianLinearDegenerate def;
  def.Init(aXYZ0, aTri, true);
  def.aBCFlag.assign(aXYZ0.size(),0);
  def.max_itr = 1000;
  def.weight_nrm = 1;
  for(unsigned int iframe=0;iframe<30;iframe++){
    // deform target mesh
    trg.Def(0.5*sin(0.1*iframe));

    // sample the space
    for(auto& smplr: sampler.aSampler){
      smplr.InitGL(); // move the sampled image to a texture
      smplr.Start();
      dfm2::opengl::SetView(smplr);
      ::glDisable(GL_POLYGON_OFFSET_FILL ); // the depth will be jazzy without this
      ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(trg.aXYZ1,trg.aElm);
      smplr.End();
    }

    // deformation
    Project(def.aIdpNrm, aXYZ1, sampler);
    def.SetBoundaryConditionToPreconditioner();
    def.Deform(aXYZ1, aXYZ0);
    std::cout << iframe << " nitr:" << def.aConvHist.size() << " nconst:" << def.aIdpNrm.size() << " np:" << aXYZ0.size()/3 << std::endl;

    // ----------------------------
    // drawing functions
    viewer.DrawBegin_oldGL();
    Draw(trg,aXYZ1,aTri);
    sampler.Draw();
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

int main()
{
  dfm2::opengl::CRender2Tex_DrawOldGL_BOX sampler;
  sampler.Initialize(128, 128, 128, 0.015);
  
  // ---------------------------------------
  CCapsuleRigged trg;
  trg.Init();
 
  // ---------------------------------------
  dfm2::glfw::CViewer3 viewer(2.0);
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);
    
  // ---------------------
  // making depth surface
  
  std::vector<double> aXYZ1;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3D_Sphere(aXYZ1, aTri, 0.7, 32, 32);
  const std::vector<double> aXYZ0 = aXYZ1; // initial src mesh vertex

  // ---------------------
  for(unsigned int itr=0;itr<3;++itr){
    aXYZ1 = aXYZ0;
    LaplacianLinear(
        aXYZ1,trg,sampler,viewer,
        aXYZ0,aTri);
    aXYZ1 = aXYZ0;
    LaplacianDegenerate(
        aXYZ1,trg,sampler,viewer,
        aXYZ0,aTri);
  }
}


