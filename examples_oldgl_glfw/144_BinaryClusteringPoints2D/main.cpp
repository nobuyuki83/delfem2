/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cstdlib>
#include <set>
#include <random>
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/color.h"
#include "delfem2/clusterpoints.h"
#include "delfem2/dtri2_v2dtri.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"

namespace dfm2 = delfem2;

// ----------------------------

void MeshTri2D_Square(
    std::vector<double>& aXY,
    std::vector<unsigned int>& aTri,
    double elen)
{
  namespace dfm2 = delfem2;
  //
  std::vector< std::vector<double> > aaXY0;
  {
    std::vector<double> aXY0;
    aXY0.push_back(0); aXY0.push_back(0);
    aXY0.push_back(1); aXY0.push_back(0);
    aXY0.push_back(1); aXY0.push_back(1);
    aXY0.push_back(0); aXY0.push_back(1);
    aaXY0.push_back(aXY0);
  }

  std::vector<dfm2::CDynPntSur> aPo2D;
  std::vector<dfm2::CVec2d> aVec2;
  std::vector<dfm2::CDynTri> aETri;
  { // initial tesselation
    std::vector<int> loopIP_ind, loopIP;
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY0);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if( elen > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     elen );
    }
    Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                   loopIP_ind,loopIP);
  }
  AssertDTri(aETri);
  AssertMeshDTri(aPo2D,aETri);
  CheckTri(aPo2D,aETri,aVec2);

  // filling inside
  if( elen > 1.0e-10 ){
    dfm2::CInputTriangulation_Uniform param(1.0);
    std::vector<int> aFlgPnt(aPo2D.size());
    std::vector<unsigned int> aFlgTri(aETri.size(),0);
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                  aVec2.size(), 0, elen, param);
  }
  AssertDTri(aETri);
  AssertMeshDTri(aPo2D,aETri);
  CheckTri(aPo2D,aETri,aVec2);

  // export to stl array
  dfm2::MeshTri2D_Export(aXY,aTri, aVec2, aETri);
}

// -----------------------------

int main(int argc,char* argv[])
{
  class CClusterData
  {
  public:
    std::vector<double> aXY; // center position of the cluster
    std::vector<double> aArea; // area of the cluster
    std::vector<unsigned int> psup_ind, psup; // connectivity of the cluster
    // below: data for visualization
    std::vector<float> aColor; // color of the cluster
    std::vector<unsigned int> map0c; // index of cluster for mesh points.
  };

  std::vector<CClusterData> aPointData;

  std::vector<unsigned int> aTri0;
  aPointData.resize(1);
  {
    {
      CClusterData &pd0 = aPointData[0];
      MeshTri2D_Square(pd0.aXY,aTri0,0.02);
    }
    const unsigned int np0 = aPointData[0].aXY.size() / 2;
    { // make area
      CClusterData &pd0 = aPointData[0];
      pd0.aArea.resize(np0);
      dfm2::MassPoint_Tri2D(pd0.aArea.data(),
                            1.0,
                            pd0.aXY.data(), np0,
                            aTri0.data(), aTri0.size() / 3);
    }
    { // make psup
      CClusterData &pd0 = aPointData[0];
      dfm2::JArray_PSuP_MeshElem(pd0.psup_ind, pd0.psup,
                                 aTri0.data(), aTri0.size() / 3, 3,
                                 np0);
    }
    { // initialize map
      aPointData[0].map0c.resize(np0);
      for(unsigned int ip=0;ip<np0;++ip){
        aPointData[0].map0c[ip] = ip;
      }
    }
  }

  for(unsigned int itr=0;itr<8;++itr) {
    aPointData.resize(aPointData.size() + 1);
    const CClusterData& pd0 = aPointData[itr];
    CClusterData& pd1 = aPointData[itr + 1];
    std::vector<unsigned int> map01(pd0.aXY.size()/2);
    unsigned int np1 = dfm2::BinaryClustering_Points2d(
        map01.data(),
        pd0.aXY.size()/2, pd0.aArea.data(), pd0.psup_ind.data(), pd0.psup.data());
    {
      const unsigned int np0 = pd0.aXY.size()/2;
      pd1.aXY.assign(np1*2,0.0);
      pd1.aArea.assign(np1,0.0);
      for(unsigned int ip0=0;ip0<np0;++ip0){
        const unsigned int ip1 = map01[ip0];
        const double a0 = pd0.aArea[ip0];
        pd1.aArea[ip1] += a0;
        pd1.aXY[ip1*2+0] += a0*pd0.aXY[ip0*2+0];
        pd1.aXY[ip1*2+1] += a0*pd0.aXY[ip0*2+1];
      }
      for(unsigned int ip1=0;ip1<np1;++ip1){
        const double a1inv = 1.0/pd1.aArea[ip1];
        pd1.aXY[ip1*2+0] *= a1inv;
        pd1.aXY[ip1*2+1] *= a1inv;
      }
    }
    dfm2::Clustering_Psup(
        pd1.psup_ind,pd1.psup,
        pd1.aXY.size()/2,
        pd0.aXY.size()/2,map01.data(),pd0.psup_ind.data(),pd0.psup.data());
    unsigned int np0 = aPointData[0].aXY.size()/2;
    pd1.map0c.resize(np0,UINT_MAX);
    for(unsigned int ip=0;ip<np0;++ip){
      unsigned int ic0 = pd0.map0c[ip];
      assert( ic0 < map01.size() );
      pd1.map0c[ip] = map01[ic0];
      assert( pd1.map0c[ip] < pd1.aXY.size()/2 );
    }
  }

  for( auto& pd : aPointData ){
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<double> dist(0, 1.0);
    const unsigned int np = pd.aXY.size()/2;
    pd.aColor.resize(np*3);
    for(unsigned int ip=0;ip<np;++ip) {
      float *pc = pd.aColor.data() + ip * 3;
      dfm2::GetRGB_HSV(pc[0], pc[1], pc[2],
                       dist(eng), 1.0, 1.0);
    }
  }

  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.nav.camera.trans[0] = -0.5;
  viewer.nav.camera.trans[1] = -0.5;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 0.7;
  while (!glfwWindowShouldClose(viewer.window) )
  {
    for(const auto& pd: aPointData) {
      const CClusterData &dp0 = aPointData[0];
      for(unsigned int itr=0;itr<30;++itr) {
        viewer.DrawBegin_oldGL();
        ::glBegin(GL_TRIANGLES);
        for (unsigned int it = 0; it < aTri0.size() / 3; ++it) {
          const unsigned int i0 = aTri0[it * 3 + 0];
          const unsigned int i1 = aTri0[it * 3 + 1];
          const unsigned int i2 = aTri0[it * 3 + 2];
          const unsigned int ic0 = pd.map0c[i0];  assert( ic0 < pd.aColor.size()/3 );
          const unsigned int ic1 = pd.map0c[i1];  assert( ic1 < pd.aColor.size()/3 );
          const unsigned int ic2 = pd.map0c[i2];  assert( ic2 < pd.aColor.size()/3 );
          ::glColor3fv(pd.aColor.data() + ic0 * 3);
          ::glVertex2dv(dp0.aXY.data() + i0 * 2);
          ::glColor3fv(pd.aColor.data() + ic1 * 3);
          ::glVertex2dv(dp0.aXY.data() + i1 * 2);
          ::glColor3fv(pd.aColor.data() + ic2 * 3);
          ::glVertex2dv(dp0.aXY.data() + i2 * 2);
        }
        ::glEnd();
        ::glColor3d(0,0,0);
        //
        /*
        dfm2::opengl::DrawMeshTri2D_Edge(aPointData[0].aXY.data(),
                                         aPointData[0].aXY.size()/2,
                                         aTri0.data(),
                                         aTri0.size()/3);
                                         */
        //
        ::glColor3d(0, 0, 0);
        ::glTranslated(0,0,+0.01);
        dfm2::opengl::DrawPoints2d_Psup(
            pd.aXY.size()/2, pd.aXY.data(),
            pd.psup_ind.data(), pd.psup.data());
        ::glTranslated(0,0,-0.01);
        viewer.SwapBuffers();
        glfwPollEvents();
      }
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
