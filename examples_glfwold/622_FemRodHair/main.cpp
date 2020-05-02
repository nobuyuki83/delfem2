/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femrod.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"
// --------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/funcs_glold.h"

namespace dfm2 = delfem2;

// -------------------------------------

void myGlutDisplay
(const std::vector<dfm2::CVec3d>& aP,
 const std::vector<dfm2::CVec3d>& aS,
 std::vector<unsigned int>& aIP_HairRoot)
{
  const unsigned int nhair = aIP_HairRoot.size()-1;
  for(unsigned int ihair=0;ihair<nhair;++ihair){
    const unsigned int ips = aIP_HairRoot[ihair];
    const unsigned int ipe = aIP_HairRoot[ihair+1];
    assert( aP.size() == aS.size() );
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1,0,0);
    ::glPointSize(3);
    ::glBegin(GL_POINTS);
    for(unsigned int ip=ips;ip<ipe;++ip){
      ::glVertex3d(aP[ip].x(), aP[ip].y(), aP[ip].z());
    }
    ::glEnd();
    // ------------
    ::glColor3d(0,0,0);
    ::glLineWidth(3);
    ::glBegin(GL_LINES);
    unsigned int ns = ipe-ips-1;
    for(unsigned int is=0;is<ns;++is){
      const unsigned int ip0 = ips+is+0; assert( ip0 < aP.size() );
      const unsigned int ip1 = ips+is+1; assert( ip1 < aP.size() );
      ::glVertex3d(aP[ip0].x(), aP[ip0].y(), aP[ip0].z());
      ::glVertex3d(aP[ip1].x(), aP[ip1].y(), aP[ip1].z());
    }
    ::glEnd();
    // --------------
    ::glBegin(GL_LINES);
    for(unsigned int is=0;is<ns;++is){
      const unsigned int ip0 = ips+is+0; assert( ip0 < aP.size() );
      const unsigned int ip1 = ips+is+1; assert( ip1 < aP.size() );
      dfm2::CVec3d p01 = 0.5*(aP[ip0]+aP[ip1]);
      double l01 = (aP[ip0]-aP[ip1]).Length();
      dfm2::opengl::myGlVertex(p01);
      dfm2::opengl::myGlVertex(p01+(l01*0.5)*aS[is]);
    }
    ::glEnd();
  }
}

class CHairShape{
public:
  unsigned int np;
  double pitch;
  double rad0;
  double dangle;
  double p0[3];
};

void MakeProblemSetting_Spiral
(std::vector<dfm2::CVec3d>& aP0,
 std::vector<dfm2::CVec3d>& aS0,
 std::vector<unsigned int>& aIP_HairRoot,
 const std::vector<CHairShape>& aHairShape)
{
  aIP_HairRoot.assign(1,0);
  aP0.clear();
  aS0.clear();
  for(unsigned int ihair=0;ihair<aHairShape.size();++ihair){
    const unsigned int np = aHairShape[ihair].np;
    const double pitch = aHairShape[ihair].pitch;
    const double dangle = aHairShape[ihair].dangle;
    const double rad0 = aHairShape[ihair].rad0;
    const double* p0 = aHairShape[ihair].p0;
    for(unsigned int ip=0;ip<np;++ip){
      dfm2::CVec3d p = dfm2::CVec3d(p0[0]+ip*pitch, p0[1]+rad0*cos(dangle*ip), p0[2]+rad0*sin(dangle*ip));
      aP0.push_back(p);
    };
    const unsigned int np0 = aIP_HairRoot[ihair];
    for(unsigned int is=0;is<np-1;++is){
      const dfm2::CVec3d v = (aP0[np0+is+1] - aP0[np0+is+0]).Normalize();
      dfm2::CVec3d s(1.3, 1.5, 1.7);
      s = (s-(s*v)*v).Normalize();
      aS0.push_back(s);
    }
    aS0.push_back(dfm2::CVec3d(1,0,0));
    aIP_HairRoot.push_back(aP0.size());
  }
}


void ParallelTransport
(std::vector<dfm2::CVec3d>& aP0,
 std::vector<dfm2::CVec3d>& aS0,
 const std::vector<unsigned int>& aIP_HairRoot)
{
  assert( aP0.size() == aS0.size() );
  assert( !aIP_HairRoot.empty() && aIP_HairRoot[0] == 0 );
  assert( aP0.size() == aIP_HairRoot[aIP_HairRoot.size()-1]);
  for(unsigned int ih=0;ih<aIP_HairRoot.size()-1;++ih){
    const unsigned int ip_r = aIP_HairRoot[ih];
    const unsigned int np = aIP_HairRoot[ih+1]-ip_r;
    for(unsigned int ir=0;ir<np-2;++ir){
      const unsigned int ip0 = ip_r+ir+0;
      const unsigned int ip1 = ip_r+ir+1;
      const unsigned int ip2 = ip_r+ir+2;
      const unsigned int is0 = ip_r+ir+0;
      const unsigned int is1 = ip_r+ir+1;
      const dfm2::CMat3d CMat3 = dfm2::Mat3_MinimumRotation(aP0[ip1]-aP0[ip0], aP0[ip2]-aP0[ip1]);
      dfm2::CVec3d s1 = CMat3*aS0[is0] + aS0[is1];
      const dfm2::CVec3d v = (aP0[ip2] - aP0[ip1]).Normalize();
      aS0[is1] = (s1-(s1*v)*v).Normalize();
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
  std::random_device rd;
  std::mt19937 reng(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  // ------
  std::vector<dfm2::CVec3d> aP0, aS0;
  std::vector<unsigned int> aIP_HairRoot;
  std::vector<int> aBCFlag;
  dfm2::CMatrixSparse<double> mats;
  std::vector<dfm2::CVec3d> aS, aP;
  // -------
  int iframe = 0;
  while (true)
  {
    if( iframe % 70 == 0 ){
      {
        std::vector<CHairShape> aHairShape;
        double rad0 = dist(reng);
        double dangle = dist(reng);
        for(int ihair=0;ihair<10;++ihair){
          CHairShape hs;
          hs.np = 30;
          hs.pitch = 0.1;
          hs.rad0 = rad0;
          hs.dangle = dangle;
          hs.p0[0] = -1;
          hs.p0[1] = (dist(reng)-0.5)*2.0;
          hs.p0[2] = (dist(reng)-0.5)*2.0;
          aHairShape.push_back(hs);
        }
        MakeProblemSetting_Spiral(aP0,aS0,aIP_HairRoot,
                                  aHairShape); // dangle
      }
      for(int itr=0;itr<10;++itr){
        ParallelTransport(aP0,aS0,aIP_HairRoot);
      }
      {
        aBCFlag.assign(aP0.size()*4, 0);
        for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
          assert( aIP_HairRoot[ihair+1] > aIP_HairRoot[ihair] );
          unsigned int ips0 = aIP_HairRoot[ihair]+0;
          unsigned int ips1 = aIP_HairRoot[ihair]+1;
          unsigned int ipe1 = aIP_HairRoot[ihair+1]-1;
          aBCFlag[ips0*4+0] = 1; aBCFlag[ips0*4+1] = 1; aBCFlag[ips0*4+2] = 1;
          aBCFlag[ips1*4+0] = 1; aBCFlag[ips1*4+1] = 1; aBCFlag[ips1*4+2] = 1;
          aBCFlag[ips0*4+3] = 1;
          aBCFlag[ipe1*4+3] = 1;
        }
      }
      {
        const unsigned int np = aP0.size();
        std::vector<unsigned int> psup_ind, psup;
        {
          const unsigned int nr = np-2;
          std::vector<unsigned int> aElemRod(nr*3);
          for(unsigned int ir=0;ir<nr;++ir){
            aElemRod[ir*3+0] = ir+0;
            aElemRod[ir*3+1] = ir+1;
            aElemRod[ir*3+2] = ir+2;
          }
          dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                                     aElemRod.data(), aElemRod.size()/3, 3, np);
        }
        dfm2::JArray_Sort(psup_ind, psup);
        mats.Initialize(np, 4, true);
        mats.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
      }
      // -----------------
      aP = aP0;
      aS = aS0;
      assert( aS.size() == aS0.size() );
      assert( aP.size() == aP0.size() );
      assert( aS.size() == aP.size() );
      // apply random deviation
      for(unsigned int ip=0;ip<aP.size();++ip){
        aP[ip] = aP0[ip];
        auto rnd = dfm2::CVec3d::Random()*0.03;
        if( aBCFlag[ip*4+0] == 0 ){ aP[ip].p[0] += rnd.x(); }
        if( aBCFlag[ip*4+1] == 0 ){ aP[ip].p[1] += rnd.y(); }
        if( aBCFlag[ip*4+2] == 0 ){ aP[ip].p[2] += rnd.z(); }
        if( aBCFlag[ip*4+3] == 0 ){
          assert( ip != aP.size()-1 );
          aS[ip] += dfm2::CVec3d::Random()*0.03;
        }
      }
      for(unsigned int is=0;is<aP.size()-1;++is){
        assert( is < aS.size() );
        const unsigned int ip0 = is+0;
        const unsigned int ip1 = is+1;
        const dfm2::CVec3d& p0 = aP[ip0];
        const dfm2::CVec3d& p1 = aP[ip1];
        const dfm2::CVec3d e01 = (p1-p0).Normalize();
        aS[is] -= (aS[is]*e01)*e01;
        aS[is].SetNormalizedVector();
      }
    }
    { // DispRotSeparate
      Solve_DispRotCombined(aP, aS, mats,
                            aP0, aS0, aBCFlag, aIP_HairRoot);
    }
    iframe = iframe+1;
    // -------------
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aP,aS,aIP_HairRoot);
    viewer.DrawEnd_oldGL();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
