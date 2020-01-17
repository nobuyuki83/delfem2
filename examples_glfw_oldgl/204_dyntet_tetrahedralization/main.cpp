/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include "delfem2/vec3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/dtet_v3.h"

// ------------

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"

// --------------------------------------------

void Inactivate
(int it,
 std::vector<CETet>& aSTet)
{
  for(int ift=0;ift<4;++ift){
    int jt0 = aSTet[it].s[ift];
    if( jt0 == -1 ){
    }
    else{
      int jft0 = tetRel[ aSTet[it].f[ift] ][ift];
      assert( aSTet[jt0].s[jft0] == it );
      aSTet[jt0].s[jft0] = -1;
      aSTet[it].s[ift] = -1;
    }
  }
  aSTet[it].v[0] = -1;
  aSTet[it].v[1] = -1;
  aSTet[it].v[2] = -1;
  aSTet[it].v[3] = -1;
}

//　-------------------------------------

std::vector<CEPo3D> aPo3D;
std::vector<CETet> aSTet;

int imode_display = 0;
bool is_animation = true;

int it_wrong = -1;

// --------------------------------------


static void myGlVertex3d(const CVector3& v)
{
  ::glVertex3d(v.x, v.y, v.z);
}


void myGlutDisplay()
{

  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  ::glEnable(GL_CULL_FACE);
  ::glCullFace(GL_BACK);
  
  
  if (imode_display==1){
    ::glColor3d(0, 0, 0);
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    myGlVertex3d(aPo3D[631].p);
    myGlVertex3d(aPo3D[732].p);
    myGlVertex3d(aPo3D[964].p);
    ::glEnd();
  }
  
//  else if(imode_display ==1){
  {
    ::glColor3d(0, 0, 0);
    ::glPointSize(1);
    ::glBegin(GL_POINTS);
    for (auto & ip : aPo3D){
      glVertex3d(ip.p.x, ip.p.y, ip.p.z);
    }
    ::glEnd();
    ///
    ::glBegin(GL_LINES);
    for (auto & it : aSTet){
      int ip0 = it.v[0];
      if( ip0 == -1 ) continue;
      int ip1 = it.v[1];
      int ip2 = it.v[2];
      int ip3 = it.v[3];
      myGlVertex3d(aPo3D[ip0].p);
      myGlVertex3d(aPo3D[ip1].p);
      myGlVertex3d(aPo3D[ip0].p);
      myGlVertex3d(aPo3D[ip2].p);
      myGlVertex3d(aPo3D[ip0].p);
      myGlVertex3d(aPo3D[ip3].p);
      myGlVertex3d(aPo3D[ip1].p);
      myGlVertex3d(aPo3D[ip2].p);
      myGlVertex3d(aPo3D[ip1].p);
      myGlVertex3d(aPo3D[ip3].p);
      myGlVertex3d(aPo3D[ip2].p);
      myGlVertex3d(aPo3D[ip3].p);
    }
    ::glEnd();
  }
  
}



void Initialize()
{
  aPo3D.clear();
  aSTet.clear();
  {
    double len = 3.0;
    aPo3D.resize(4);
    aPo3D[0].p = CVector3(-len, +len, -len); aPo3D[0].e = 0; aPo3D[0].poel = 0;
    aPo3D[1].p = CVector3(+len, -len, -len); aPo3D[1].e = 0; aPo3D[1].poel = 1;
    aPo3D[2].p = CVector3(+len, +len, +len); aPo3D[2].e = 0; aPo3D[2].poel = 2;
    aPo3D[3].p = CVector3(-len, -len, +len); aPo3D[3].e = 0; aPo3D[3].poel = 3;
    aSTet.resize(1);
    aSTet[0].v[0] = 0;
    aSTet[0].v[1] = 1;
    aSTet[0].v[2] = 2;
    aSTet[0].v[3] = 3;
    aSTet[0].s[0] = -1;
    aSTet[0].s[1] = -1;
    aSTet[0].s[2] = -1;
    aSTet[0].s[3] = -1;
    aSTet[0].setCircumCenter(aPo3D);
  }

}

void AddRandomPoint()
{
  std::vector<int> tmp_buffer;
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<> dist(-1.0, +1.0);
  for(int itr=0;itr<100;itr++){
    double x0 = dist(mt);
    double y0 = dist(mt);
    double z0 = dist(mt);
    int ip_ins = (int)aPo3D.size();
    aPo3D.resize(ip_ins+1);
    aPo3D[ip_ins].p = CVector3(x0,y0,z0);
    aPo3D[ip_ins].e = -1;
    aPo3D[ip_ins].poel = -1;
    int itet_ins = -1;
    {
      for (std::size_t it = 0; it<aSTet.size(); ++it){
        int j0 = aSTet[it].v[0];
        int j1 = aSTet[it].v[1];
        int j2 = aSTet[it].v[2];
        int j3 = aSTet[it].v[3];
        if( j0 == -1 ) continue; // floating tet
        double v0 = TetVolume(ip_ins, j1, j2, j3, aPo3D);
        double v1 = TetVolume(j0, ip_ins, j2, j3, aPo3D);
        double v2 = TetVolume(j0, j1, ip_ins, j3, aPo3D);
        double v3 = TetVolume(j0, j1, j2, ip_ins, aPo3D);
        //    double v4 = TetVolume(j0, j1, j2, j3, aPo3D);
        if (v0>0&&v1>0&&v2>0&&v3>0){ itet_ins = it; break; }
      }
    }
    if (itet_ins==-1){ return; }
    AddPointTetDelaunay(ip_ins,itet_ins, aPo3D, aSTet, tmp_buffer);
#ifndef NDEBUG
    CheckTet(aSTet, aPo3D);
#endif
  }
  std::cout << aPo3D.size() << std::endl;
}


void myGlutKeyboard(unsigned char Key, int x, int y)
{
    {
      std::vector<double> aXYZ;
      std::vector<unsigned int> aTri;
      delfem2::Read_Obj("models/bunny2k.obj",
                        aXYZ, aTri);
      delfem2::Normalize_Points3D(aXYZ);
      delfem2::Scale_PointsX(aXYZ,
                             2.0);
      // -------------------------
      Initialize();
      std::vector<int> tmp_buffer;
      for(std::size_t ixyz=0;ixyz<aXYZ.size()/3;ixyz++){
        double x0 = aXYZ[ixyz*3+0];
        double y0 = aXYZ[ixyz*3+1];
        double z0 = aXYZ[ixyz*3+2];
        int ip_ins = (int)aPo3D.size();
        aPo3D.resize(ip_ins+1);
        aPo3D[ip_ins].p = CVector3(x0,y0,z0);
        aPo3D[ip_ins].e = -1;
        aPo3D[ip_ins].poel = -1;
        int itet_ins = -1;
        { // find tet
          for (std::size_t it = 0; it<aSTet.size(); ++it){
            int j0 = aSTet[it].v[0];
            int j1 = aSTet[it].v[1];
            int j2 = aSTet[it].v[2];
            int j3 = aSTet[it].v[3];
            if( j0 == -1 ) continue; // floating tet
            double v0 = TetVolume(ip_ins, j1, j2, j3, aPo3D);
            double v1 = TetVolume(j0, ip_ins, j2, j3, aPo3D);
            double v2 = TetVolume(j0, j1, ip_ins, j3, aPo3D);
            double v3 = TetVolume(j0, j1, j2, ip_ins, aPo3D);
            //    double v4 = TetVolume(j0, j1, j2, j3, aPo3D);
            if (v0>-0.0000001&&v1>-0.00000001&&v2>-0.000000001&&v3>-0.00000001){ itet_ins = it; break; }
          }
        }
        if (itet_ins==-1){ continue; }
        AddPointTetDelaunay(ip_ins,itet_ins, aPo3D, aSTet, tmp_buffer);
#ifndef NDEBUG
        CheckTet(aSTet, aPo3D);
#endif
      }
      { // remove SuperTet verteces
        for(std::size_t it=0;it<aSTet.size();++it){
          int i0 = aSTet[it].v[0];
          int i1 = aSTet[it].v[1];
          int i2 = aSTet[it].v[2];
          int i3 = aSTet[it].v[3];
          if( i0 < 4 || i1 < 4 || i2 < 4 || i3 < 4 ){
            Inactivate(it,aSTet);
          }
        }
        aPo3D[0].e = -1;
        aPo3D[1].e = -1;
        aPo3D[2].e = -1;
        aPo3D[3].e = -1;
        for(std::size_t it=0;it<aSTet.size();++it){
          for(int ift=0;ift<4;++ift){
            int iv = aSTet[it].v[ift];
            if( iv == -1 ) continue;
            aPo3D[iv].e = it;
            aPo3D[iv].poel = ift;
          }
        }
      }
      { // edge recovery
        std::vector<unsigned int> psup_ind, psup;
        JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                                    aTri.data(),aTri.size()/3,3,(int)aXYZ.size()/3);
        std::vector<unsigned int> edge_ind, edge;
        JArrayEdgeUnidir_PointSurPoint(edge_ind, edge,
                                psup_ind, psup);
//        CJaggedArray edge;
//        edge.SetEdgeOfElem(aTri, (int)aTri.size()/3, 3, (int)aXYZ.size()/3, false);
        for(int ixyz=0;ixyz<(int)aXYZ.size()/3;++ixyz){
          int ip0 = ixyz+4;
          ElemAroundPoint elarpo;
          {
            int itet0 = aPo3D[ip0].e;
            if( itet0 == -1 ){
              std::cout << ip0 << " " << aPo3D[ip0].e << std::endl;
              continue;
            }
            MakeElemAroundPoint(elarpo, itet0, aPo3D[ip0].poel, aSTet);
          }
          for(int ie=edge_ind[ixyz];ie<edge_ind[ixyz+1];++ie){
            int jxyz = edge[ie];
            int jp0 = jxyz+3;
          }
        }
      }
    }
}


int main(int argc, char* argv[])
{
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();

  viewer.nav.camera.view_height = 2.5;
  
  Initialize();
  
  while(!glfwWindowShouldClose(viewer.window)){
    {
      static int iframe = 0;
      if( iframe == 0 ){
        Initialize();
      }
      else if( iframe % 100 == 0 ){
        AddRandomPoint();
      }
      iframe = (iframe+1)%500;
    }
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
