/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <vector>
#include "delfem2/dtri_v2.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/glold_v23dtricad.h"

namespace dfm2 = delfem2;

// --------------------------------------------

std::vector<dfm2::CDynPntSur> aPo2D;
std::vector<CVector2> aVec2;
std::vector<dfm2::CDynTri> aETri;
std::vector<int> loopIP_ind, loopIP;

// --------------------------------------------

void Refine(double px, double py)
{
  dfm2::CCmdRefineMesh aCmd;
  RefinementPlan_EdgeLongerThan_InsideCircle(aCmd,
                                             0.05, px, py, 0.1,
                                             aPo2D, aVec2, aETri);
  RefineMesh(aPo2D, aETri, aVec2, aCmd);
}

void Coarse(double px, double py)
{
  for(int ip=(int)aPo2D.size()-1;ip>=0;--ip){
    if( aPo2D[ip].e == -1 ){ continue; }
    if( Distance(aVec2[ip],CVector2(px,py)) > 0.1 ){ continue; }
    std::vector< std::pair<int,int> > aTriSuP;
    GetTriArrayAroundPoint(aTriSuP,
                           ip,aPo2D,aETri);
    std::vector<int> aPSuP(aTriSuP.size());
    const int npsup = aPSuP.size();
    for(int iit=0;iit<npsup;++iit){
      const int itri0 = aTriSuP[iit].first;
      const int inotri0 = aTriSuP[iit].second;
      assert( ip == aETri[itri0].v[inotri0] );
      aPSuP[iit] = aETri[itri0].v[ (inotri0+2)%3 ];
    }
    std::map<double,int> mapDistTri;
    for(int iit=0;iit<npsup;++iit){
      int ip1 = aPSuP[iit];
      double d01 = Distance(aVec2[ip],aVec2[ip1]);
      double min_area = 0.0;
      for(int jjt=0;jjt<npsup-2;++jjt){
        const int ip2 = aPSuP[(iit+jjt+1)%npsup];
        const int ip3 = aPSuP[(iit+jjt+2)%npsup];
        double area = TriArea(aVec2[ip1],aVec2[ip2],aVec2[ip3]);
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
        const int itri0 = aTriSuP[iit0].first;
        const int inotri0 = aTriSuP[iit0].second;
        Collapse_ElemEdge(itri0, (inotri0+1)%3, aPo2D, aETri);
        const int ip1 = aETri[itri0].v[ (inotri0+2)%3 ];
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
  /////
  const double elen = 0.11;
  {
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if( elen > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     elen );
    }
  }
  ////
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind,loopIP);
  if( elen > 1.0e-10 ){
    dfm2::CInputTriangulation_Uniform param(1.0);
    std::vector<int> aFlgPnt(aPo2D.size()), aFlgTri(aETri.size());
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
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
      ::glVertex3d(aVec2[ip].x(), aVec2[ip].y(), 0.1);
    }
    ::glEnd();
  }
}

int main(int argc,char* argv[])
{
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
    
  GenMesh();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      static int iframe = 0;
      if( iframe == 0 ){
        GenMesh();
      }
      if( iframe % 100 == 0 ){
        double px = (rand()/(RAND_MAX+1.0))-0.5;
        double py = (rand()/(RAND_MAX+1.0))-0.5;
        if( iframe < 2000 ){
          Refine(px, py);
        }
        else{
          Coarse(px, py);
        }
      }
      iframe = (iframe + 1)%4000;
    }
    // ----------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
