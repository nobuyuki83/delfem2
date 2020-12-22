/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/femrod.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/mats.h"
#include "delfem2/mshuni.h"
#include "delfem2/jagarray.h"
#include <GLFW/glfw3.h>


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
  for(const auto & p : aP){
    ::glVertex3d(p.x(), p.y(), p.z());
  }
  ::glEnd();
  // ------------
  ::glColor3d(0,0,0);
  ::glLineWidth(3);
  ::glBegin(GL_LINES);
  for(unsigned int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    unsigned int i0 = aElemSeg[iseg*2+0]; assert( i0 < aP.size() );
    unsigned int i1 = aElemSeg[iseg*2+1]; assert( i1 < aP.size() );
    ::glVertex3d(aP[i0].x(), aP[i0].y(), aP[i0].z());
    ::glVertex3d(aP[i1].x(), aP[i1].y(), aP[i1].z());
  }
  ::glEnd();
  // --------------
  ::glBegin(GL_LINES);
  for(unsigned int iseg=0;iseg<aElemSeg.size()/2;++iseg){
    unsigned int i0 = aElemSeg[iseg*2+0]; assert( i0 < aP.size() );
    unsigned int i1 = aElemSeg[iseg*2+1]; assert( i1 < aP.size() );
    dfm2::CVec3d p01 = 0.5*(aP[i0]+aP[i1]);
    double l01 = (aP[i0]-aP[i1]).Length();
    dfm2::opengl::myGlVertex(p01);
    dfm2::opengl::myGlVertex(p01+(l01*0.5)*aS[iseg]);
  }
  ::glEnd();
}

void MakeProblemSetting_Spiral
(std::vector<dfm2::CVec3d>& aP0,
 std::vector<dfm2::CVec3d>& aS0,
 std::vector<unsigned int>& aElemSeg,
 std::vector<unsigned int>& aElemRod,
 unsigned int np,
 double pitch,
 double rad0,
 double dangle)
{
  aP0.resize(np);
  for(unsigned int ip=0;ip<np;++ip){
    aP0[ip] = dfm2::CVec3d(-1.0+ip*pitch, rad0*cos(dangle*ip), rad0*sin(dangle*ip));
  };
  // -------------------------
  // below: par segment data
  const unsigned int ns = np-1;
  aElemSeg.resize(ns*2);
  for(unsigned int is=0;is<ns;++is){
    aElemSeg[is*2+0] = is+0;
    aElemSeg[is*2+1] = is+1;
  }
  { // initial director vector
    aS0.resize(ns,dfm2::CVec3d(1,0,0));
    for(unsigned int is=0;is<ns;++is){
      unsigned int ip0 = aElemSeg[is*2+0];
      unsigned int ip1 = aElemSeg[is*2+1];
      const dfm2::CVec3d v = (aP0[ip1] - aP0[ip0]).Normalize();
      aS0[is] = (aS0[is]-(aS0[is]*v)*v).Normalize();
    }
  }
  // --------------------------
  // below: par rod element data
  const unsigned int nr = ns-1;
  aElemRod.resize(nr*5);
  for(unsigned int ir=0;ir<nr;++ir){
    aElemRod[ir*5+0] = ir+0;
    aElemRod[ir*5+1] = ir+1;
    aElemRod[ir*5+2] = ir+2;
    aElemRod[ir*5+3] = np+ir+0;
    aElemRod[ir*5+4] = np+ir+1;
  };
  // smoothing
  for(int itr=0;itr<10;++itr){
    for(unsigned int ir=0;ir<nr;++ir){
      const unsigned int ip0 = aElemRod[ir*5+0];
      const unsigned int ip1 = aElemRod[ir*5+1];
      const unsigned int ip2 = aElemRod[ir*5+2];
      const unsigned int is0 = aElemRod[ir*5+3]-np; assert( is0 < ns );
      const unsigned int is1 = aElemRod[ir*5+4]-np; assert( is1 < ns );
      const dfm2::CMat3d CMat3 = dfm2::Mat3_MinimumRotation(aP0[ip1]-aP0[ip0], aP0[ip2]-aP0[ip1]);
      dfm2::CVec3d s1 = CMat3*aS0[is0] + aS0[is1];
      const dfm2::CVec3d v = (aP0[ip2] - aP0[ip1]).Normalize();
      aS0[is1] = (s1-(s1*v)*v).Normalize();
    }
  }
}


int main(int argc,char* argv[])
{
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  // -----
  std::random_device rd;
  std::mt19937 reng(rd());
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  std::uniform_real_distribution<double> dist03(0.0, 3.0);
  // ------
  while (true)
  {
    std::vector<dfm2::CVec3d> aP0, aS0;
    std::vector<unsigned int> aElemSeg, aElemRod;
    MakeProblemSetting_Spiral(aP0, aS0,
                              aElemSeg, aElemRod,
                              30,
                              0.1, // np
                              dist01(reng), // rad0
                              dist01(reng)); // dangle
    std::vector<int> aBCFlag;
    {
      const unsigned int np = aP0.size();
      const unsigned int ns = aS0.size();
      const unsigned int nNode = np+ns;
      aBCFlag.assign(nNode*3, 0);
      {
        aBCFlag[0*3+0] = 1; aBCFlag[0*3+1] = 1; aBCFlag[0*3+2] = 1;
        aBCFlag[1*3+0] = 1; aBCFlag[1*3+1] = 1; aBCFlag[1*3+2] = 1;
        aBCFlag[(np+0)*3+0] = 1; //
        for(unsigned int is=0;is<ns;++is){
          aBCFlag[(np+is)*3+1] = 1; // fix the unused dof
          aBCFlag[(np+is)*3+2] = 1; // fix the unused dof
        }
      }
    }
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
    std::vector<dfm2::CVec3d> aS = aS0, aP = aP0;
    assert( aS.size() == aElemSeg.size()/2);
    // apply random deviation
    for(unsigned int ip=0;ip<aP.size();++ip){
      aP[ip] = aP0[ip];
      auto rnd = dfm2::CVec3d::Random(dist03,reng);
      if( aBCFlag[ip*3+0] == 0 ){ aP[ip].p[0] += rnd.x(); }
      if( aBCFlag[ip*3+1] == 0 ){ aP[ip].p[1] += rnd.y(); }
      if( aBCFlag[ip*3+2] == 0 ){ aP[ip].p[2] += rnd.z(); }
    }
    const unsigned int ns = aS.size();
    for(unsigned int is=0;is<ns;++is){
      aS[is] = aS0[is];
      auto rnd = dfm2::CVec3d::Random(dist03,reng);
      const unsigned int np = aP.size();
      if( aBCFlag[(np+is)*3+0] == 0 ){ aS[is].p[0] += rnd.x(); }
      if( aBCFlag[(np+is)*3+1] == 0 ){ aS[is].p[1] += rnd.y(); }
      if( aBCFlag[(np+is)*3+2] == 0 ){ aS[is].p[2] += rnd.z(); }
    }
    for(unsigned int iseg=0;iseg<aElemSeg.size()/2;++iseg){
      const unsigned int i0 = aElemSeg[iseg*2+0];
      const unsigned int i1 = aElemSeg[iseg*2+1];
      const dfm2::CVec3d& p0 = aP[i0];
      const dfm2::CVec3d& p1 = aP[i1];
      const dfm2::CVec3d e01 = (p1-p0).Normalize();
      assert( iseg < aS.size() );
      aS[iseg] -= (aS[iseg]*e01)*e01;
      aS[iseg].SetNormalizedVector();
    }
    const double stiff_stretch = dist01(reng)+1.;
    const double stiff_bendtwist[3] = {
        dist01(reng)+1.,
        dist01(reng)+1.,
        dist01(reng)+1. };
    for(int iframe=0;iframe<70;++iframe) {
      Solve_DispRotSeparate(
          aP, aS, mats,
          stiff_stretch, stiff_bendtwist,
          aP0, aS0, aElemSeg, aElemRod, aBCFlag);
      viewer.DrawBegin_oldGL();
      myGlutDisplay(aP, aS, aElemSeg);
      viewer.SwapBuffers();
      glfwPollEvents();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
   }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
