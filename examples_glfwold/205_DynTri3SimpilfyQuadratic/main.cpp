/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/dtri3_v3dtri.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -----------------------------

void myGlutDisplay
 (const std::vector<dfm2::CDynPntSur>& aPo,
  const std::vector<dfm2::CDynTri>& aTri,
  const std::vector<dfm2::CVec3d>& aVec3)
{
  GLboolean is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glEnable(GL_LIGHTING);
  GLboolean is_texture  = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_TEXTURE_2D);    
  {
    //    float gray[4] = {0.3,0.3,0.3,1};
    float gray[4] = {0.9f,0.9f,0.9f,1.f};    
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    float shine[4] = {0,0,0,0};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
    ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 127.0);
    //    ::glColor3d(1,1,1);
  }
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(const auto & tri : aTri){
    const unsigned int i1 = tri.v[0];
    const unsigned int i2 = tri.v[1];
    const unsigned int i3 = tri.v[2];
    dfm2::opengl::myGlVertex(aVec3[i1]);
    dfm2::opengl::myGlVertex(aVec3[i2]);
    dfm2::opengl::myGlVertex(aVec3[i3]);
  }
  ::glEnd();        
  
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(const auto & tri : aTri){
    const unsigned int i1 = tri.v[0];
    const unsigned int i2 = tri.v[1];
    const unsigned int i3 = tri.v[2];
    dfm2::opengl::myGlVertex(aVec3[i1]);
    dfm2::opengl::myGlVertex(aVec3[i2]);
    dfm2::opengl::myGlVertex(aVec3[i2]);
    dfm2::opengl::myGlVertex(aVec3[i3]);
    dfm2::opengl::myGlVertex(aVec3[i3]);
    dfm2::opengl::myGlVertex(aVec3[i1]);
  }
  ::glEnd();
  // ----------------
  if( is_lighting ){ ::glEnable( GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }
  if(  is_texture  ){ ::glEnable(GL_TEXTURE_2D); }
}

double MinimizeQuad(double* pos,
                    const double* SymMat4_Q)
{
  const dfm2::CMat3d A(SymMat4_Q[0], SymMat4_Q[1], SymMat4_Q[2],
                       SymMat4_Q[1], SymMat4_Q[4], SymMat4_Q[5],
                       SymMat4_Q[2], SymMat4_Q[5], SymMat4_Q[7]);
  const dfm2::CVec3d v(SymMat4_Q[3], SymMat4_Q[6], SymMat4_Q[8]);
  const dfm2::CVec3d p = -A.Inverse()*v;
  p.CopyTo(pos);
  return SymMat4_Q[9] + p*v;
//  std::cout << ip << " ###  " << aVec3[ip] << "  ###  " << p - aVec3[ip] << " ### " << err << std::endl;
}


void QuadMetric
 (std::vector<double>& aSymMat4,
  std::vector<dfm2::CDynPntSur>& aDP,
  std::vector<dfm2::CDynTri>& aDTri,
  std::vector<dfm2::CVec3d>& aVec3)
{
  const unsigned int np = aDP.size();
  aSymMat4.assign(np*10, 0.0);
  for(unsigned int ip=0;ip<np;++ip){
    std::vector< std::pair<unsigned int,unsigned int> > aTriSurPo;
    dfm2::GetTriArrayAroundPoint(aTriSurPo,
                                 ip,aDP,aDTri);
    for(unsigned iit=0;iit<aTriSurPo.size();++iit){
      unsigned int it0 = aTriSurPo[iit].first;
      unsigned int ino0 = aTriSurPo[iit].second;
      assert( aDTri[it0].v[ino0] == ip );
      const dfm2::CVec3d n0 = dfm2::UnitNormal_DTri3(it0, aDTri, aVec3);
      const double a0 = n0.x();
      const double b0 = n0.y();
      const double c0 = n0.z();
      const double d0 = -n0*aVec3[ip];
      aSymMat4[ip*10+0] += a0*a0;
      aSymMat4[ip*10+1] += a0*b0;
      aSymMat4[ip*10+2] += a0*c0;
      aSymMat4[ip*10+3] += a0*d0;
      aSymMat4[ip*10+4] += b0*b0;
      aSymMat4[ip*10+5] += b0*c0;
      aSymMat4[ip*10+6] += b0*d0;
      aSymMat4[ip*10+7] += c0*c0;
      aSymMat4[ip*10+8] += c0*d0;
      aSymMat4[ip*10+9] += d0*d0;
    }
    double pos[3];
    double err = MinimizeQuad(pos, aSymMat4.data()+ip*10);
    std::cout << ip << " ### " << err << " ### " << dfm2::CVec3d(pos)-aVec3[ip] << std::endl;
  }
  
}

int main(int argc,char* argv[])
{
  std::vector<dfm2::CDynPntSur> aDP;
  std::vector<dfm2::CDynTri> aDTri;
  std::vector<dfm2::CVec3d> aVec3;
  {
    std::vector<double> aXYZ0;
    std::vector<unsigned int> aTri0;
    delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/arm_16k.ply",
                      aXYZ0,aTri0);
    dfm2::Normalize_Points3(aXYZ0,2.0);
    const unsigned int np = aXYZ0.size()/3;
    aDP.resize(np);
    aVec3.resize(np);
    for(unsigned int ipo=0;ipo<aDP.size();ipo++){
      aVec3[ipo].p[0] = aXYZ0[ipo*3+0];
      aVec3[ipo].p[1] = aXYZ0[ipo*3+1];
      aVec3[ipo].p[2] = aXYZ0[ipo*3+2];
    }
    InitializeMesh(aDP, aDTri,
                   aTri0.data(),aTri0.size()/3,aVec3.size());
#if !defined(NDEBUG)
    AssertDTri(aDTri);
    AssertMeshDTri(aDP, aDTri);
    AssertMeshDTri2(aDP, aDTri, aVec3);
#endif
  }
  std::vector<double> aSymMat4;
  QuadMetric(aSymMat4,
             aDP,aDTri,aVec3);
  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.5;
  while (!glfwWindowShouldClose(viewer.window))
  {
    for(unsigned int i=0;i<10;i++){
      auto itri0 = (unsigned int)((rand()/(RAND_MAX+1.0))*aDTri.size());
      assert( itri0 < aDTri.size() );
      const int ip_del = aDTri[itri0].v[2];
      bool res = CollapseElemEdge(itri0, 0, aDP, aDTri);
#if !defined(NDEBUG)
      if( res ){
        assert( aDP[ip_del].e == -1 );
      }
      AssertMeshDTri(aDP, aDTri);
#endif
      if( aDTri.size() <= 100 ) break;
    }
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aDP,aDTri,aVec3);
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
