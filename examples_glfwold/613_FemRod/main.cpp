/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cstdlib>
#include <vector>
#include <set>
#include <random>
#include "delfem2/mat3.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/objfunc_v23.h"
//
#include "delfem2/v23m34q.h"

// --------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/glold_funcs.h"

namespace dfm2 = delfem2;

// -------------------------------------

void myGlutDisplay
(const std::vector<dfm2::CVec3d>& aP,
 const std::vector<dfm2::CVec3d>& aS,
 const std::vector<unsigned int>& aElemSeg)
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for(int ip=0;ip<aP.size();++ip){
    ::glVertex3d(aP[ip].x(), aP[ip].y(), aP[ip].z());
  }
  ::glEnd();
  // ------------
  ::glColor3d(0,0,0);
  ::glLineWidth(3);
  ::glBegin(GL_LINES);
  for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    unsigned int i0 = aElemSeg[iseg*2+0]; assert( i0 < aP.size() );
    unsigned int i1 = aElemSeg[iseg*2+1]; assert( i1 < aP.size() );
    ::glVertex3d(aP[i0].x(), aP[i0].y(), aP[i0].z());
    ::glVertex3d(aP[i1].x(), aP[i1].y(), aP[i1].z());
  }
  ::glEnd();
  // --------------
  ::glBegin(GL_LINES);
  for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    unsigned int i0 = aElemSeg[iseg*2+0]; assert( i0 < aP.size() );
    unsigned int i1 = aElemSeg[iseg*2+1]; assert( i1 < aP.size() );
    dfm2::CVec3d p01 = 0.5*(aP[i0]+aP[i1]);
    double l01 = (aP[i0]-aP[i1]).Length();
    dfm2::opengl::myGlVertex(p01);
    dfm2::opengl::myGlVertex(p01+(l01*0.5)*aS[iseg]);
  }
  ::glEnd();
}

void Solve
(std::vector<dfm2::CVec3d>& aP,
 std::vector<dfm2::CVec3d>& aS,
 dfm2::CMatrixSparse<double>& mats,
 const std::vector<dfm2::CVec3d>& aP0,
 const std::vector<unsigned int>& aElemSeg,
 const std::vector<unsigned int>& aElemRod,
 const std::vector<int>& aBCFlag)
{
  const unsigned int nNode = aBCFlag.size()/3;
  mats.SetZero();
  std::vector<double> vec_r;
  vec_r.assign(nNode*3, 0.0);
  std::vector<int> tmp_buffer;
  double W = 0;
  for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    const unsigned int i0 = aElemSeg[iseg*2+0];
    const unsigned int i1 = aElemSeg[iseg*2+1];
    const unsigned int* aINoel = aElemSeg.data()+iseg*2;
    const double L0 = (aP0[i0]-aP0[i1]).Length();
    const dfm2::CVec3d aPE[2] = { aP[i0], aP[i1] };
    // --------------
    dfm2::CVec3d dW_dP[2];
    dfm2::CMat3d ddW_ddP[2][2];
    W += WdWddW_SquareLengthLineseg3D(dW_dP, ddW_ddP,
                                      aPE, L0);
    {
      double eM[2*2*3*3];
      for(int in=0;in<2;++in){
        for(int jn=0;jn<2;++jn){
          ddW_ddP[in][jn].CopyTo(eM+(in*2+jn)*9);
        }
      }
      mats.Mearge(2, aINoel, 2, aINoel, 9, eM, tmp_buffer);
    }
    {
      for (int inoel=0; inoel<2; inoel++){
        const unsigned int ip = aINoel[inoel];
        vec_r[ip*3+0] -= dW_dP[inoel].x();
        vec_r[ip*3+1] -= dW_dP[inoel].y();
        vec_r[ip*3+2] -= dW_dP[inoel].z();
      }
    }
  }
  for(int irod=0;irod<aElemRod.size()/5;++irod){
    const unsigned int* aINoel = aElemRod.data()+irod*5;
    const unsigned int nP = aP.size();
    const dfm2::CVec3d aPE[3] = { aP[aINoel[0]], aP[aINoel[1]], aP[aINoel[2]] };
    const dfm2::CVec3d aSE[2] = {
      aS[aINoel[3]-nP],
      aS[aINoel[4]-nP] };
    const double off[3] = {0,0,0};
    // ------
    dfm2::CVec3d dW_dP[3];
    double dW_dt[2];
    dfm2::CMat3d ddW_ddP[3][3];
    dfm2::CVec3d ddW_dtdP[2][3];
    double ddW_ddt[2][2];
    W +=  WdWddW_Rod(dW_dP,dW_dt,ddW_ddP,ddW_dtdP,ddW_ddt,
                     aPE,aSE,off, false);
    {
      double eM[5][5][3][3];
      for(int i=0;i<5*5*3*3;++i){ (&eM[0][0][0][0])[i] = 0.0; }
      for(int in=0;in<3;++in){
        for(int jn=0;jn<3;++jn){
//          dfm2::CMat3d t = ddW_ddP[in][jn].Trans();
//          t.CopyValueToPtr(&eM[in][jn][0][0]);
          ddW_ddP[in][jn].CopyTo(&eM[in][jn][0][0]);
        }
      }
      for(int in=0;in<2;++in){
        for(int jn=0;jn<3;++jn){
          eM[3+in][jn][0][0] = ddW_dtdP[in][jn].x();
          eM[3+in][jn][0][1] = ddW_dtdP[in][jn].y();
          eM[3+in][jn][0][2] = ddW_dtdP[in][jn].z();
        }
      }
      for(int in=0;in<3;++in){
        for(int jn=0;jn<2;++jn){
          eM[in][jn+3][0][0] = ddW_dtdP[jn][in].x();
          eM[in][jn+3][1][0] = ddW_dtdP[jn][in].y();
          eM[in][jn+3][2][0] = ddW_dtdP[jn][in].z();
        }
      }
      for(int in=0;in<2;++in){
        for(int jn=0;jn<2;++jn){
          eM[in+3][jn+3][0][0] = ddW_ddt[in][jn];
        }
      }
      mats.Mearge(5, aINoel, 5, aINoel, 9, &eM[0][0][0][0], tmp_buffer);
    }
    {
      for (int inoel=0; inoel<3; inoel++){
        const unsigned int ip = aINoel[inoel];
        vec_r[ip*3+0] -= dW_dP[inoel].x();
        vec_r[ip*3+1] -= dW_dP[inoel].y();
        vec_r[ip*3+2] -= dW_dP[inoel].z();
      }
      for (int inoel=0; inoel<2; inoel++){
        const unsigned int in0 = aINoel[3+inoel];
        vec_r[in0*3+0] -= dW_dt[inoel];
      }
    }
  }
//  std::cout << dfm2::CheckSymmetry(mats) << std::endl;
//  mats.AddDia(0.1);
  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << dfm2::CheckSymmetry(mats) << std::endl;
  mats.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_r, aBCFlag,0);
  std::vector<double> vec_x;
  vec_x.assign(nNode*3, 0.0);
 {
   auto aConvHist = Solve_CG(vec_r.data(),vec_x.data(),
   vec_r.size(), 1.0e-4, 100, mats);
   if( aConvHist.size() > 0 ){
     std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
   }
 }
  /*
  {
    auto aConvHist = dfm2::Solve_BiCGStab(vec_r,vec_x,
                                          1.0e-4, 300, mats);
    if( aConvHist.size() > 0 ){
      std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
    }
  }
   */
  //    for(int i=0;i<vec_x.size();++i){
  //      std::cout << i << " " << vec_x[i] << std::endl;
  //    }
  assert( aS.size() == aElemSeg.size()/2 );
  for(int is=0;is<aS.size();++is){
    unsigned int i0 = aElemSeg[is*2+0];
    unsigned int i1 = aElemSeg[is*2+1];
    dfm2::CVec3d V01 = aP[i1]-aP[i0];
    dfm2::CVec3d du(vec_x[i1*3+0]-vec_x[i0*3+0],
                    vec_x[i1*3+1]-vec_x[i0*3+1],
                    vec_x[i1*3+2]-vec_x[i0*3+2]);
    const unsigned int np = aP.size();
    const double dtheta = vec_x[ np*3 + is*3 ];
    dfm2::CVec3d frm[3];
    RodFrameTrans(frm,
                  aS[is],V01,du,dtheta);
    aS[is] = frm[0];
  }
  for(int ip=0;ip<aP.size();++ip){
    aP[ip].p[0] += vec_x[ip*3+0];
    aP[ip].p[1] += vec_x[ip*3+1];
    aP[ip].p[2] += vec_x[ip*3+2];
  }
  for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    const unsigned int i0 = aElemSeg[iseg*2+0];
    const unsigned int i1 = aElemSeg[iseg*2+1];
    const dfm2::CVec3d& p0 = aP[i0];
    const dfm2::CVec3d& p1 = aP[i1];
    const dfm2::CVec3d e01 = (p1-p0).Normalize();
    aS[iseg] -= (aS[iseg]*e01)*e01;
    aS[iseg].SetNormalizedVector();
  }
}


void MakeProblemSetting
(std::vector<dfm2::CVec3d>& aP0,
 std::vector<unsigned int>& aElemSeg,
 std::vector<unsigned int>& aElemRod,
 std::vector<int>& aBCFlag, // if value this is not 0, it is fixed boundary condition
 unsigned int np)
{
  double elen = 2.0/np;
  aP0.resize(np);
  for(int ip=0;ip<np;++ip){
    aP0[ip] = dfm2::CVec3d(-1.0+ip*elen, 0.0, 0.0);
  };
  unsigned int ns = np-1;
  aElemSeg.resize(ns*2);
  for(int is=0;is<ns;++is){
    aElemSeg[is*2+0] = is+0;
    aElemSeg[is*2+1] = is+1;
  }
  unsigned int nr = ns-1;
  aElemRod.resize(nr*5);
  for(int ir=0;ir<nr;++ir){
    aElemRod[ir*5+0] = ir+0;
    aElemRod[ir*5+1] = ir+1;
    aElemRod[ir*5+2] = ir+2;
    aElemRod[ir*5+3] = np+ir+0;
    aElemRod[ir*5+4] = np+ir+1;
  };
  const unsigned int nNode = np+ns;
  aBCFlag.assign(nNode*3, 0);
  {
    unsigned int np = aP0.size();
    aBCFlag[0*3+0] = 1; aBCFlag[0*3+1] = 1; aBCFlag[0*3+2] = 1;
    aBCFlag[1*3+0] = 1; aBCFlag[1*3+1] = 1; aBCFlag[1*3+2] = 1;
    aBCFlag[(np+0)*3+0] = 1;
    for(int is=0;is<ns;++is){
      aBCFlag[(np+is)*3+1] = 1;
      aBCFlag[(np+is)*3+2] = 1;
    }
  }
}

int main(int argc,char* argv[])
{
  // -----
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.5;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  // -----
  std::vector<dfm2::CVec3d> aP0;
  std::vector<unsigned int> aElemSeg, aElemRod;
  std::vector<int> aBCFlag;
  MakeProblemSetting(aP0,aElemSeg,aElemRod, aBCFlag,
                     30);
  dfm2::CMatrixSparse<double> mats;
  {
    unsigned int nNode = aElemSeg.size()/2 + aP0.size();
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aElemRod.data(), aElemRod.size()/5, 5, nNode);
    dfm2::JArray_Sort(psup_ind, psup);
    mats.Initialize(nNode, 3, true);
    mats.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
  }
  // -----------------
  std::vector<dfm2::CVec3d> aP = aP0;
  std::vector<dfm2::CVec3d> aS(aElemSeg.size()/2);
  assert( aP.size() == aP0.size() );
  assert( aS.size() == aElemSeg.size()/2);
  // ------
  int iframe = 0;
  while (true)
  {
    if( iframe % 200 == 0 ){
      for(int ip=0;ip<aP.size();++ip){
        aP[ip] = aP0[ip];
        auto rnd = dfm2::CVec3d::Random()*3.0;
        if( aBCFlag[ip*3+0] == 0 ){ aP[ip].p[0] += rnd.x(); }
        if( aBCFlag[ip*3+1] == 0 ){ aP[ip].p[1] += rnd.y(); }
        if( aBCFlag[ip*3+2] == 0 ){ aP[ip].p[2] += rnd.z(); }
      }
      unsigned int ns = aS.size();
      for(int is=0;is<ns;++is){
        aS[is] = dfm2::CVec3d(0,1,0);
        auto rnd = dfm2::CVec3d::Random()*3.0;
        const unsigned int np = aP.size();
        if( aBCFlag[(np+is)*3+0] == 0 ){ aS[is].p[0] += rnd.x(); }
        if( aBCFlag[(np+is)*3+1] == 0 ){ aS[is].p[1] += rnd.y(); }
        if( aBCFlag[(np+is)*3+2] == 0 ){ aS[is].p[2] += rnd.z(); }
      }
      for(int iseg=0;iseg<aElemSeg.size()/2;++iseg){
        const unsigned int i0 = aElemSeg[iseg*2+0];
        const unsigned int i1 = aElemSeg[iseg*2+1];
        const dfm2::CVec3d& p0 = aP[i0];
        const dfm2::CVec3d& p1 = aP[i1];
        const dfm2::CVec3d e01 = (p1-p0).Normalize();
        assert( iseg < aS.size() );
        aS[iseg] -= (aS[iseg]*e01)*e01;
        aS[iseg].SetNormalizedVector();
      }
    }
    {
      Solve(aP, aS, mats,
            aP0, aElemSeg, aElemRod, aBCFlag);
    }
    iframe = iframe+1;
    // -------------
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aP,aS,aElemSeg);
    viewer.DrawEnd_oldGL();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
