/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cmath>
#include <vector>
#include "delfem2/dtri2_v2dtri.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/cad2dtriv2_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// --------------------------------------

double AreaCGCurve(const std::vector<double>& aCV, double cg[2])
{
  const unsigned int nCV = aCV.size()/2;
  double area = 0;
  cg[0] = 0;
  cg[1] = 0;
  for(unsigned int idiv=0;idiv<nCV;idiv++){
    unsigned int ipo0 = idiv;
    unsigned int ipo1 = idiv+1;
    if( idiv == nCV-1 ){ ipo1 = 0; }
    const double p0[2] = { aCV[ipo0*2+0], aCV[ipo0*2+1] };
    const double p1[2] = { aCV[ipo1*2+0], aCV[ipo1*2+1] };
    const double p2[2] = { 0, 0 };
    double a0 = dfm2::Area_Tri2(p0, p1, p2);
    double cg0[2] = { (p0[0]+p1[0]+p2[0])/3.0, (p0[1]+p1[1]+p2[1])/3.0 };
    cg[0] += cg0[0]*a0;
    cg[1] += cg0[1]*a0;
    area += a0;
  }
  cg[0] /= area;
  cg[1] /= area;
  return area;
}

void MakeRandomCV(unsigned int nCV, std::vector<double>& aCV)
{
  aCV.clear();
  for(unsigned int icv=0;icv<nCV;icv++){
    /*
     {
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     }
     */
    {
      double tht = icv*3.1415*2.0/nCV;
      double r = (double)rand()/(RAND_MAX+1.0);
      double px = r*sin(tht);
      double py = r*cos(tht);
      aCV.push_back(px);
      aCV.push_back(py);
    }
  }
  
  {
    double cnt[2];
    double area = AreaCGCurve(aCV,cnt);
    if( area < 0 ){
      std::vector<double> aCV0;
      for(unsigned int idiv=0;idiv<nCV;idiv++){
        unsigned int idiv0 = nCV-1-idiv;
        aCV0.push_back( aCV[idiv0*2+0] );
        aCV0.push_back( aCV[idiv0*2+1] );
      }
      aCV = aCV0;
      area = -area;
    }
    if( area < 1.0e-5 ){
      aCV.clear();
      return;
    }
  }
}

void MakeCurveSpline(const std::vector<double>& aCV, std::vector<double>& aVecCurve)
{
  aVecCurve.resize(0);
  const unsigned int nCV = aCV.size()/2;
  unsigned int ndiv = 5;
  for(unsigned int icv=0;icv<nCV;icv++){
    unsigned int icv0=(icv+0)%nCV;
    unsigned int icv1=(icv+1)%nCV;
    unsigned int icv2=(icv+2)%nCV;
    const double p0[2] = { aCV[icv0*2+0], aCV[icv0*2+1] };
    const double p1[2] = { aCV[icv1*2+0], aCV[icv1*2+1] };
    const double p2[2] = { aCV[icv2*2+0], aCV[icv2*2+1] };
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      const double t = 1.0-(double)idiv/ndiv;
      const double w[3] = {0.5*t*t, -t*t + t + 0.5, 0.5*(1-t)*(1-t) };
      const double px = w[0]*p0[0] + w[1]*p1[0] + w[2]*p2[0];
      const double py = w[0]*p0[1] + w[1]*p1[1] + w[2]*p2[1];
      aVecCurve.push_back(px);
      aVecCurve.push_back(py);
    }
  }
}

void myGlVertex2D(const std::vector<double>& vec, unsigned int i)
{
  ::glVertex3d(vec[i*2],vec[i*2+1],+0.5);
}


void drawCurve
(const std::vector<double>& vec,
 const std::vector<double>& aVecCurve0)
{
  if( aVecCurve0.size() < 2 ){ return; }
  ::glBegin(GL_LINES);
  const unsigned int nvec = vec.size()/2;
  for(unsigned int ivec=0;ivec<nvec;ivec++){
    unsigned int jvec = ivec+1; if( jvec >= nvec ){ jvec -= nvec; }
    myGlVertex2D(vec,ivec);
    myGlVertex2D(vec,jvec);
  }
  ::glEnd();
  
  ::glBegin(GL_POINTS);
  for(unsigned int ivec=0;ivec<nvec;ivec++){
    myGlVertex2D(vec,ivec);
  }
  ::glEnd();
}

void drawMesh
(std::vector<int>& aTri,
 std::vector<double>& aXY)
{
  const unsigned int ntri = aTri.size()/3;
  const unsigned int nxys = aXY.size()/2;
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  //  double mag = 20;
  for(unsigned int itri=0;itri<ntri;itri++){
    const int ino0 = aTri[itri*3+0];
    const int ino1 = aTri[itri*3+1];
    const int ino2 = aTri[itri*3+2];
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
  }
  ::glEnd();
  ////////////////
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itri=0;itri<ntri;itri++){
    const int ino0 = aTri[itri*3+0];
    const int ino1 = aTri[itri*3+1];
    const int ino2 = aTri[itri*3+2];
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
  }
  ::glEnd();
  //
  ::glPointSize(2);
  ::glColor3d(0,0,0);
  ::glBegin(GL_POINTS);
  for(unsigned int ino=0;ino<nxys;ino++){
    ::glVertex2d( aXY[ino*2+0], aXY[ino*2+1] );
  }
  ::glEnd();
  
}



// --------------------------------

std::vector<dfm2::CDynPntSur> aPo2D;
std::vector<dfm2::CVec2d> aVec2;
std::vector<dfm2::CDynTri> aETri;
std::vector<int> loopIP_ind, loopIP;

int idp_nearest = -1;

int press_button = -1;
double mov_begin_x, mov_begin_y;
bool is_animation = true;
double mag = 1.0;

// --------------------------

void GenMesh(){
  std::vector<double> aCV0; MakeRandomCV(8,aCV0); // current cv
  std::vector<double> aVecCurve0;  MakeCurveSpline(aCV0,aVecCurve0); // current curve
  ////
  std::vector< std::vector<double> > aaXY;
  aaXY.push_back( aVecCurve0 );
  /////
  const double elen = 0.03;
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
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind,loopIP);
  if( elen > 1.0e-10 ){
    dfm2::CInputTriangulation_Uniform param(1.0);
    std::vector<int> aFlgPnt(aPo2D.size()), aFlgTri(aETri.size());
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                  aVec2.size(), 0, elen, param);
  }
}

// --------------------------------------

void myGlutDisplay(void)
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
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      static int iframe = 0;
      if( iframe == 0 ){
        GenMesh();
      }
      iframe = (iframe + 1)%300;
    }
    
    viewer.DrawBegin_oldGL();
    
    myGlutDisplay();
    
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
