/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <vector>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/cad2dtriv2.h"

namespace dfm2 = delfem2;

// --------------------------------------------

std::vector<dfm2::CDynPntSur> aPo2D;
std::vector<dfm2::CVec2d> aVec2;
std::vector<dfm2::CDynTri> aETri;
std::vector<int> loopIP_ind, loopIP;

// --------------------------------------------

void Refine(double px, double py)
{
  dfm2::CCmdRefineMesh aCmd;
  RefinementPlan_EdgeLongerThan_InsideCircle(
      aCmd,
      0.05, px, py, 0.1,
      aPo2D, aVec2, aETri);
  RefineMesh(aPo2D, aETri, aVec2, aCmd);
}

void Coarse(double px, double py)
{
  for(auto ip=static_cast<unsigned int>(aPo2D.size());ip-->0;){ // iterate ip-1 to 0
    if( aPo2D[ip].e == UINT_MAX ){ continue; }
    if( Distance(aVec2[ip],dfm2::CVec2d(px,py)) > 0.1 ){ continue; }
    std::vector< std::pair<unsigned int,unsigned int> > aTriSuP;
    GetTriArrayAroundPoint(aTriSuP,
                           ip,aPo2D,aETri);
    std::vector<int> aPSuP(aTriSuP.size());
    const auto npsup = static_cast<unsigned int>(aPSuP.size());
    for(unsigned int iit=0;iit<npsup;++iit){
      const unsigned int itri0 = aTriSuP[iit].first;
      const unsigned int inotri0 = aTriSuP[iit].second;
      assert( ip == aETri[itri0].v[inotri0] );
      aPSuP[iit] = aETri[itri0].v[ (inotri0+2)%3 ];
    }
    std::map<double,int> mapDistTri;
    for(unsigned int iit=0;iit<npsup;++iit){
      int ip1 = aPSuP[iit];
      double d01 = Distance(aVec2[ip],aVec2[ip1]);
      double min_area = 0.0;
      for(unsigned int jjt=0;jjt<npsup-2;++jjt){
        const int ip2 = aPSuP[(iit+jjt+1)%npsup];
        const int ip3 = aPSuP[(iit+jjt+2)%npsup];
        double area = Area_Tri(aVec2[ip1],aVec2[ip2],aVec2[ip3]);
        if( jjt == 0 || area < min_area ){ min_area = area; }
      }
      if( min_area > 1.0e-10 ){
        mapDistTri.insert( std::make_pair(d01,iit) );
      }
    }
    if( !mapDistTri.empty() ){
      const int iit0 = mapDistTri.begin()->second;
      assert( iit0>=0 && iit0 < (int)aPSuP.size() );
      double dist0 = mapDistTri.begin()->first;
      if( dist0 < 0.05 ){
        const unsigned int itri0 = aTriSuP[iit0].first;
        const unsigned int inotri0 = aTriSuP[iit0].second;
        CollapseEdge_MeshDTri(itri0, (inotri0+1)%3, aPo2D, aETri);
        const unsigned int ip1 = aETri[itri0].v[ (inotri0+2)%3 ];
        DelaunayAroundPoint(ip1, aPo2D, aETri, aVec2);
      }
    }
  }
}


void GenMesh()
{
  std::vector< std::vector<double> > aaXY;
  aaXY.resize(1);
  {
    double xys[8] = {-0.5,-0.5, +0.5,-0.5, +0.5,+0.5, -0.5,+0.5};
    aaXY[0].assign(xys,xys+8);
  }
  // --------------------------------
  constexpr double elen = 0.11;
  {
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if constexpr( elen > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     elen );
    }
  }
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind,loopIP);
  if constexpr( elen > 1.0e-10 ){
    dfm2::CInputTriangulation_Uniform param(1.0);
    std::vector<int> aFlgPnt(aPo2D.size());
    std::vector<unsigned int> aFlgTri(aETri.size(),0);
    MeshingInside(
        aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
        aVec2.size(), 0, elen, param);
  }
}

void myGlutDisplay()
{
  ::glPointSize(5);
  ::glLineWidth(1);
  ::glPointSize(5);
  ::glColor3d(1,1,0);
  delfem2::opengl::DrawMeshDynTri_Edge(aETri, aVec2);
  ::glColor3d(0.8, 0.8, 0.8);
  delfem2::opengl::DrawMeshDynTri_FaceNorm(aETri, aVec2);
  ::glLineWidth(3);
  ::glColor3d(0,0,0);
  for(int iloop=0;iloop<(int)loopIP_ind.size()-1;iloop++){
    ::glBegin(GL_LINE_LOOP);
    for(int iip=loopIP_ind[iloop];iip<loopIP_ind[iloop+1];iip++){
      const int ip = loopIP[iip];
      ::glVertex3d(aVec2[ip].x, aVec2[ip].y, 0.1);
    }
    ::glEnd();
  }
}

int main()
{
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();

  std::mt19937 random_engine(std::random_device{}());
  std::uniform_real_distribution<double> dist_m05p05(-0.5, +0.5);
    
  GenMesh();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      static int iframe = 0;
      if( iframe == 0 ){
        GenMesh();
      }
      else{
        double px = dist_m05p05(random_engine);
        double py = dist_m05p05(random_engine);
        if( iframe < 40 ){
          Refine(px, py);
        }
        else{
          Coarse(px, py);
        }
      }
      iframe = (iframe + 1)%40;
    }
    // ----------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
