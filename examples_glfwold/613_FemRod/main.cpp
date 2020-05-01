/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femrod.h"
#include "delfem2/objf_geo3.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"
#include "delfem2/geo3_v23m34q.h"


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
 const std::vector<unsigned int>& aElemSeg)
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for(unsigned int ip=0;ip<aP.size();++ip){
    ::glVertex3d(aP[ip].x(), aP[ip].y(), aP[ip].z());
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
  std::vector<dfm2::CVec3d> aP0, aS0, aDarboux0;
  std::vector<unsigned int> aElemSeg, aElemRod;
  std::vector<int> aBCFlag;
  dfm2::CMatrixSparse<double> mats;
  std::vector<dfm2::CVec3d> aS, aP;
  // -------
  int iframe = 0;
  while (true)
  {
    if( iframe % 70 == 0 ){
      MakeProblemSetting_Spiral(aP0,aS0,aDarboux0,
                                aElemSeg,aElemRod, aBCFlag,
                                30,
                                0.1,
                                dist(reng),
                                dist(reng));
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
      aP = aP0;
      aS = aS0;
      assert( aP.size() == aP0.size() );
      assert( aS.size() == aElemSeg.size()/2);
      // apply random deviation
      for(unsigned int ip=0;ip<aP.size();++ip){
        aP[ip] = aP0[ip];
        auto rnd = dfm2::CVec3d::Random()*3.0;
        if( aBCFlag[ip*3+0] == 0 ){ aP[ip].p[0] += rnd.x(); }
        if( aBCFlag[ip*3+1] == 0 ){ aP[ip].p[1] += rnd.y(); }
        if( aBCFlag[ip*3+2] == 0 ){ aP[ip].p[2] += rnd.z(); }
      }
      const unsigned int ns = aS.size();
      for(unsigned int is=0;is<ns;++is){
        aS[is] = aS0[is];
        auto rnd = dfm2::CVec3d::Random()*3.0;
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
    }
    {
      Solve_DispRotSeparate(aP, aS, mats,
            aP0, aDarboux0, aElemSeg, aElemRod, aBCFlag);
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
