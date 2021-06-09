/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/dtri3_v3dtri.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <climits>


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


void QuadErrorMetric_MeshDTri3
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
      const double a0 = n0.x;
      const double b0 = n0.y;
      const double c0 = n0.z;
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

class CollapseSchedule {
public:
  CollapseSchedule(unsigned int iv1,
                   unsigned int iv2,
                   const double pos[3]): iv1(iv1), iv2(iv2)
  {
    p[0] = pos[0];
    p[1] = pos[1];
    p[2] = pos[2];
  }
public:
  unsigned int iv1, iv2;
  double p[3];
};

void RemoveOnePoint
 (std::vector<dfm2::CDynPntSur>& aDP,
  std::vector<dfm2::CDynTri>& aDTri,
  std::vector<dfm2::CVec3d>& aVec3,
  std::map<double, CollapseSchedule>& cost2edge,
  std::map<std::pair<unsigned int, unsigned int>,double>& edge2cost,
  std::vector<double> aSymMat4)
{
  unsigned int itri0;
  unsigned int ied0 = UINT_MAX;
  dfm2::CVec3d pos;
  while(!cost2edge.empty()){
    unsigned int iv1=UINT_MAX, iv2=UINT_MAX;
    double err0;
    {
      auto itr_c2e = cost2edge.begin();
      iv1 = itr_c2e->second.iv1;
      iv2 = itr_c2e->second.iv2;
      pos = dfm2::CVec3d(itr_c2e->second.p);
      assert( iv1 < iv2 && iv1 < aDP.size() && iv2 < aDP.size() );
      cost2edge.erase(itr_c2e);
      err0 = itr_c2e->first;
    }
    if( aDP[iv1].e == UINT_MAX ){ continue; } // deleted vtx
    if( aDP[iv2].e == UINT_MAX ){ continue; } // deleted vtx
    {
      auto v12 = std::make_pair(iv1,iv2);
      auto itr_e2c = edge2cost.find(v12);
      double err1 = itr_e2c->second;
      if( fabs(err1-err0) > 1.0e-20 ){ continue; } // this edge is updated
      edge2cost.erase(itr_e2c);
    }
    unsigned int ino0, ino1;
    bool res = dfm2::FindEdge_LookAroundPoint(itri0, ino0, ino1,
                                              iv1, iv2, aDP, aDTri);
    if( !res ){ continue; }
    assert( aDTri[itri0].v[ino0] == iv1 );
    assert( aDTri[itri0].v[ino1] == iv2 );
    ied0 = 3-ino0-ino1;
    break;
  }
  if( ied0 == UINT_MAX ) return;
  // --------
  const unsigned int ip_sty = aDTri[itri0].v[(ied0+1)%3];
  const unsigned int ip_del = aDTri[itri0].v[(ied0+2)%3];
  bool res = CollapseEdge_MeshDTri(itri0, ied0, aDP, aDTri);
  if( !res ){ return; }
#if !defined(NDEBUG)
  if( res ){ assert( aDP[ip_del].e == UINT_MAX ); }
  AssertDTri(aDTri);
  AssertMeshDTri(aDP, aDTri);
#endif
  aVec3[ip_sty] = pos;
  {
    const double* Q1 = aSymMat4.data()+ip_sty*10;
    const double* Q2 = aSymMat4.data()+ip_del*10;
    double Q12[10];
    for(int i=0;i<10;++i){ Q12[i] = Q1[i] + Q2[i]; }
    for(int i=0;i<10;++i){ aSymMat4[ip_sty*10+i] = Q12[i]; }
  }
  std::vector<unsigned int> aIP;
  dfm2::FindPointAroundPoint(aIP,
                             ip_sty,aDP,aDTri);
  for(unsigned int iip=0;iip<aIP.size();++iip){
    unsigned int jp0 = aIP[iip];
    unsigned int iv1, iv2;
    if( ip_sty < jp0 ){
      iv1 = ip_sty;
      iv2 = jp0;
    }
    else{
      iv1 = jp0;
      iv2 = ip_sty;
    }
    assert( iv1 < iv2 );
    const double* Q1 = aSymMat4.data()+iv1*10;
    const double* Q2 = aSymMat4.data()+iv2*10;
    double Q12[10];
    for(unsigned int i=0;i<10;++i){ Q12[i] = Q1[i] + Q2[i]; }
    double pos[3];
    double err = MinimizeQuad(pos, Q12);
    auto v12 = std::make_pair(iv1,iv2);
    CollapseSchedule cs(iv1,iv2,pos);
    cost2edge.insert(std::make_pair(err,cs));
    edge2cost.insert(std::make_pair(v12,err));
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
  QuadErrorMetric_MeshDTri3(aSymMat4,
             aDP,aDTri,aVec3);
  std::map<double, CollapseSchedule> cost2edge;
  std::map<std::pair<unsigned int, unsigned int>,double> edge2cost;
  {
    for(unsigned int itri=0;itri<aDTri.size();++itri){
      for(unsigned int ied=0;ied<3;++ied){
        unsigned int iv1 = aDTri[itri].v[(ied+1)%3];
        unsigned int iv2 = aDTri[itri].v[(ied+2)%3];
        if( iv1 > iv2 ){ continue; }
        const double* Q1 = aSymMat4.data()+iv1*10;
        const double* Q2 = aSymMat4.data()+iv2*10;
        double Q12[10];
        for(unsigned int i=0;i<10;++i){ Q12[i] = Q1[i] + Q2[i]; }
        double pos[3];
        double err = MinimizeQuad(pos, Q12);
        auto v12 = std::make_pair(iv1,iv2);
        CollapseSchedule cs(iv1,iv2,pos);
        cost2edge.insert(std::make_pair(err,cs));
        edge2cost.insert(std::make_pair(v12,err));
      }
    }
  }
  // -----------
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.5;
  while (!glfwWindowShouldClose(viewer.window))
  {
    RemoveOnePoint(aDP,aDTri,aVec3,
                   cost2edge,edge2cost,
                   aSymMat4);
    // --------
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aDP,aDTri,aVec3);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
