/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cassert>
#include <vector>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/vec3.h"
#include "delfem2/geo_parametric.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

// -----------------------------------

dfm2::CVec3d GetPointSurf(
  double u, double v,
  int isurf,
  std::vector<int>& aIndCP,
  std::vector<dfm2::CVec3d>& aCP)
{
  int i00 = aIndCP[isurf*16+ 0];
  int i01 = aIndCP[isurf*16+ 1];
  int i02 = aIndCP[isurf*16+ 2];
  int i03 = aIndCP[isurf*16+ 3];
  int i04 = aIndCP[isurf*16+ 4];
  int i05 = aIndCP[isurf*16+ 5];
  int i06 = aIndCP[isurf*16+ 6];
  int i07 = aIndCP[isurf*16+ 7];
  int i08 = aIndCP[isurf*16+ 8];
  int i09 = aIndCP[isurf*16+ 9];
  int i10 = aIndCP[isurf*16+10];
  int i11 = aIndCP[isurf*16+11];
  int i12 = aIndCP[isurf*16+12];
  int i13 = aIndCP[isurf*16+13];
  int i14 = aIndCP[isurf*16+14];
  int i15 = aIndCP[isurf*16+15];
  return dfm2::getPointSurfaceBezierCubic(
    u,v,
    aCP[i00], aCP[i01], aCP[i02], aCP[i03],
    aCP[i04], aCP[i05], aCP[i06], aCP[i07],
    aCP[i08], aCP[i09], aCP[i10], aCP[i11],
    aCP[i12], aCP[i13], aCP[i14], aCP[i15]);
}

void AddQuads(
  std::vector<dfm2::CVec3d>& aPQuad,
  int n,
  int isurf,
  std::vector<int>& aIndCP,
  std::vector<dfm2::CVec3d>& aCP)
{
  for (int i = 0; i<(n+1); i++){
    for (int j = 0; j<(n+1); j++){
      for(int ic=0;ic<4;++ic){
        int iu = i;
        int jv = j;
        if(ic==1){ iu++; }
        if(ic==2){ iu++; jv++; }
        if(ic==3){ jv++; }
        double u = (double)iu/n;
        double v = (double)jv/n;
        dfm2::CVec3d p = GetPointSurf(u, v, isurf, aIndCP, aCP);
        aPQuad.push_back(p);
      }
    }
  }
}


// -----------------------------

std::vector<int> aIndCP;
std::vector<dfm2::CVec3d> aCP;
std::vector<dfm2::CVec3d> aPQuad;
int n = 20;

// ------------------------------

void Random()
{
  int nCP = 16;
  aCP.resize(nCP);
  {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<> dist(-1.0, 1.0);
    for(int iCP=0;iCP<16;iCP++) {
      aCP[iCP].p[0] = dist(mt);
      aCP[iCP].p[1] = dist(mt);
      aCP[iCP].p[2] = dist(mt);
    }
  }
  aIndCP.resize(16);
  for(int i=0;i<16;++i){ aIndCP[i] = i; }
  aCP[ 0] += dfm2::CVec3d(-2,-2,0);
  aCP[ 1] += dfm2::CVec3d(-2,-1,0);
  aCP[ 2] += dfm2::CVec3d(-2,+1,0);
  aCP[ 3] += dfm2::CVec3d(-2,+2,0);
  //
  aCP[ 4] += dfm2::CVec3d(-1,-2,0);
  aCP[ 5] += dfm2::CVec3d(-1,-1,0);
  aCP[ 6] += dfm2::CVec3d(-1,+1,0);
  aCP[ 7] += dfm2::CVec3d(-1,+2,0);
  //
  aCP[ 8] += dfm2::CVec3d(+1,-2,0);
  aCP[ 9] += dfm2::CVec3d(+1,-1,0);
  aCP[10] += dfm2::CVec3d(+1,+1,0);
  aCP[11] += dfm2::CVec3d(+1,+2,0);
  //
  aCP[12] += dfm2::CVec3d(+2,-2,0);
  aCP[13] += dfm2::CVec3d(+2,-1,0);
  aCP[14] += dfm2::CVec3d(+2,+1,0);
  aCP[15] += dfm2::CVec3d(+2,+2,0);
  //
  aPQuad.clear();
  AddQuads(aPQuad,n,0,aIndCP,aCP);
}


static void myGlVertex3d(const dfm2::CVec3d& v)
{
  ::glVertex3d(v.x,v.y,v.z);
}

static void myGlVertex3d(int i, const std::vector<dfm2::CVec3d>& aV)
{
  const dfm2::CVec3d& v = aV[i];
  ::glVertex3d(v.x,v.y,v.z);
}


void myGlutDisplay()
{
  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  ::glEnable(GL_CULL_FACE);
  ::glCullFace(GL_BACK);
  {
    ::glColor3d(0,0,0);
    ::glPointSize(3);
    ::glBegin(GL_POINTS);
    for(const auto & icp : aCP){
      myGlVertex3d(icp);
    }
    ::glEnd();
  }
  
  {
    int nq = (int)aPQuad.size()/4;
    ::glColor4d(1,1,1,1.0);
    ::glBegin(GL_QUADS);
    for(int iq=0;iq<nq;iq++){
      myGlVertex3d(iq*4+0,aPQuad);
      myGlVertex3d(iq*4+1,aPQuad);
      myGlVertex3d(iq*4+2,aPQuad);
      myGlVertex3d(iq*4+3,aPQuad);
    }
    ::glEnd();
    ////
    ::glColor4d(0,0,0,1.0);
    for(int iq=0;iq<nq;iq++){
      ::glBegin(GL_LINE_STRIP);
      myGlVertex3d(iq*4+0,aPQuad);
      myGlVertex3d(iq*4+1,aPQuad);
      myGlVertex3d(iq*4+2,aPQuad);
      myGlVertex3d(iq*4+3,aPQuad);
      myGlVertex3d(iq*4+0,aPQuad);
      ::glEnd();
    }
  }
}


int main()
{
  delfem2::glfw::CViewer3 viewer(4);
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  //
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      static int iframe = 0;
      if( iframe == 0 ){
        Random();
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
