/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <queue>
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vec3.h"

// gl related includes
#include <GLFW/glfw3.h>
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl24_tex.h"
#include "delfem2/opengl/glfw_viewer.hpp"

namespace dfm2 = delfem2;

// ------------------------------------------------------

class CNode
{
public:
  CNode(){
    itype =0;
    mesh_dist=10000000.0;
    weight=0;
  }
  CNode(const CNode& lhs){
    itype = lhs.itype;
    mesh_dist = lhs.mesh_dist;
    weight = lhs.weight;
  }
public:
  int itype; // 0:none 1:dist_in 2:fix
  double mesh_dist;
  double weight;
};
class CIndNode
{
public:
  CIndNode(unsigned int ino, double mesh_dist) :ino(ino), mesh_dist(mesh_dist){}
  bool operator < (const CIndNode& lhs) const {
    return this->mesh_dist < lhs.mesh_dist;
  }
public:
  int ino;
  double mesh_dist;
};

void GetLocalExpMap(
    double& u, double& v,
    const double p[3],
    const double lc[9],
    const double q[3])
{
  const double dist = Distance3D(p, q);
  if( dist < 1.0e-20 ){ u=0; v=0; return; }
  double pq[3] = { q[0]-p[0], q[1]-p[1], q[2]-p[2] };
  {
    const double dot = Dot3D(pq,lc);
    pq[0] -= dot*lc[0];
    pq[1] -= dot*lc[1];
    pq[2] -= dot*lc[2];
  }
  const double len = Length3D(pq);
  if( len < 1.0e-10 ){ u=0; v=0; return; }
  {
    const double c = dist/len;
    pq[0] *= c;
    pq[1] *= c;
    pq[2] *= c;
  }
  u = Dot3D(pq,lc+3);
  v = Dot3D(pq,lc+6);
}

void GetGeodesic(
    double& u_pq, double& v_pq,
    // -----------
    const double r[3], const double lcr[9],
    const double u_pr, const double v_pr,
    const double q[3])
{
  double u_rq0, v_rq0;
  GetLocalExpMap(u_rq0,v_rq0,   r,lcr,   q);
  //  const double u_rq1 = u_rq0*rot[0] + v_rq0*rot[1];
  //  const double v_rq1 = u_rq0*rot[2] + v_rq0*rot[3];
  u_pq = u_pr + u_rq0;
  v_pq = v_pr + v_rq0;
}

void MakeLocalCoord(
    double lrx[3], double lry[3],
    const double lrn[3],
    const double lcp[9] )
{
  double anrp[3];
  Cross3D(anrp, lrn, lcp);
  const double snrp = Length3D(anrp);
  if( snrp > 1.0e-5 ){
    const double cnrp = Dot3D(lcp,lrn);
    const double t = atan2(snrp,cnrp);
    anrp[0] *= 1.0/snrp;
    anrp[1] *= 1.0/snrp;
    anrp[2] *= 1.0/snrp;
    double rot[9];
    GetRotMatrix_Rodrigues3D(rot, anrp, t);
    VecMat3(lrx, lcp+3, rot);
    VecMat3(lry, lcp+6, rot);
  }
  else{
    lrx[0]=lcp[3];  lrx[1]=lcp[4];  lrx[2]=lcp[5];
    lry[0]=lcp[6];  lry[1]=lcp[7];  lry[2]=lcp[8];
  }
}

void MakeExpMap_Point(
    std::vector<double>& aTex,
    std::vector<double>& aLocCoord,
    int iker,
    const std::vector<double>& aXYZ,
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup)
{
  const unsigned int nXYZ = aXYZ.size()/3;
  std::vector<CNode> aNode;
  aNode.resize(nXYZ);
  aTex.assign(nXYZ*2,0.0);
  // -------
  double LocCoord_ker[9];
  for(unsigned int i=0;i<6;i++){ LocCoord_ker[i] = aLocCoord[iker*6+i]; }
  Cross3D(LocCoord_ker+6, LocCoord_ker+0, LocCoord_ker+3);
  // --------
  aNode[iker].itype = 1;
  aNode[iker].mesh_dist  = 0;
  aNode[iker].weight = 1;
  std::priority_queue<CIndNode> que;
  que.push(CIndNode(iker,0));
  while(!que.empty()){
    unsigned int ifix = que.top().ino;
    que.pop();
    assert( ifix >= 0 && ifix < nXYZ );
    assert( aNode[ifix].itype != 0 );
    if( aNode[ifix].itype != 1 ) continue;
    aNode[ifix].itype = 2;
    const double invw = 1.0/aNode[ifix].weight;
    aTex[ifix*2+0] *= invw;
    aTex[ifix*2+1] *= invw;
    // --------
    double LocCoord_fix[9];
    LocCoord_fix[0]=aLocCoord[ifix*6+0];
    LocCoord_fix[1]=aLocCoord[ifix*6+1];
    LocCoord_fix[2]=aLocCoord[ifix*6+2];
    MakeLocalCoord(LocCoord_fix+3,
                   LocCoord_fix+6,
                   LocCoord_fix+0,
                   LocCoord_ker);
    aLocCoord[ifix*6+3]=LocCoord_fix[3];
    aLocCoord[ifix*6+4]=LocCoord_fix[4];
    aLocCoord[ifix*6+5]=LocCoord_fix[5];
    // --------
    for(unsigned int ipsup=psup_ind[ifix];ipsup<psup_ind[ifix+1];++ipsup){
      unsigned int ino1 = psup[ipsup];
      if( aNode[ino1].itype == 2 ) continue;
      const double len = Distance3D(aXYZ.data()+ifix*3,aXYZ.data()+ino1*3);
      const double weight = 1.0/len;
      const double dist0 = aNode[ifix].mesh_dist+len;
      double u, v;
      GetGeodesic(u,v,
       aXYZ.data()+ifix*3,  LocCoord_fix,
       aTex[ifix*2+0],
       aTex[ifix*2+1],
       aXYZ.data()+ino1*3);
      if( aNode[ino1].itype == 0 || dist0 < aNode[ino1].mesh_dist ){ aNode[ino1].mesh_dist  = dist0; }
      que.push( CIndNode(ino1,-aNode[ino1].mesh_dist) );
      aNode[ino1].itype = 1;
      aNode[ino1].weight += weight;
      aTex[ino1*2+0] += u*weight;
      aTex[ino1*2+1] += v*weight;
    }
  }
}

// ---------------------------------------------------

std::vector<double> aXYZ;
std::vector<unsigned int> aTri;
std::vector<double> aTex;
std::vector<double> aLocCoord;
std::vector<unsigned int> psup_ind, psup;
int iker = 0;
bool is_animation = true;
int m_texName;
    
// ---------------------------------------------------

void SetNewProblem()
{
  const unsigned int nprob = 2;
  static unsigned int iprob = 0;
  
  {
    delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_2k.ply",
                      aXYZ, aTri);
    {
      double cx,cy,cz, wx,wy,wz;
      delfem2::CenterWidth_Points3D(cx,cy,cz,
                                    wx,wy,wz,
                                    aXYZ);
      delfem2::Translate_Points3D(aXYZ,
                         -cx,-cy,-cz);
      double wm = wx;
      wm = ( wx > wm ) ? wx : wm;
      wm = ( wy > wm ) ? wy : wm;
      wm = ( wz > wm ) ? wz : wm;
      delfem2::Scale_PointsXD(aXYZ,
                              2.0/wm);
    }
  }
  {
    std::vector<unsigned int> elsup_ind, elsup;
    dfm2::JArrayElemSurPoint_MeshElem(elsup_ind, elsup,
        aTri.data(), aTri.size()/3, 3, aXYZ.size()/3);
    dfm2::JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
        aTri.data(), elsup_ind, elsup, 3, aXYZ.size()/3);
  }
  {
    std::vector<double> aNorm(aXYZ.size());
    delfem2::Normal_MeshTri3D(aNorm.data(),
        aXYZ.data(), aXYZ.size()/3, aTri.data(),aTri.size()/3);
    const unsigned int np = aXYZ.size()/3;
    aLocCoord.resize(np*6);
    for(int ip=0;ip<np;++ip){
      double tmp[3];
      GetVertical2Vector3D(aNorm.data()+ip*3, aLocCoord.data()+ip*6+3, tmp);
      aLocCoord[ip*6+0] = aNorm[ip*3+0];
      aLocCoord[ip*6+1] = aNorm[ip*3+1];
      aLocCoord[ip*6+2] = aNorm[ip*3+2];
    }
    aTex.resize(np*2);
    for(int ip=0;ip<np;++ip){
      aTex[ip*2+0] = aXYZ[ip*3+0];
      aTex[ip*2+1] = aXYZ[ip*3+1];
    }
  }
  
  iprob++;
  if( iprob == nprob ){ iprob = 0; }
}

// -----------------------------------------------------

void myGlutDisplay()
{
  GLboolean is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glEnable(GL_LIGHTING);
  GLboolean is_texture  = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_TEXTURE_2D);
  
  { // ball
    ::glDisable(GL_TEXTURE_2D);
    {
      //    float gray[4] = {0.3,0.3,0.3,1};
      float gray[4] = {1.0f,0.0f,0.0f,1.f};
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
      float shine[4] = {0,0,0,0};
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
      ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 127.0);
      //    ::glColor3d(1,1,1);
    }
    double px = aXYZ[iker*3+0];
    double py = aXYZ[iker*3+1];
    double pz = aXYZ[iker*3+2];
    ::glDisable(GL_LIGHTING);
    ::glPushMatrix();
    ::glTranslated(+px, +py, +pz);
    ::glScaled(0.05,0.05,0.05);
    ::glColor3d(1,0,0);
    dfm2::opengl::DrawSphere(16,16);
    ::glPopMatrix();
    ////
    CVector3 ex(aLocCoord[iker*6+3],aLocCoord[iker*6+4],aLocCoord[iker*6+5]);
    ex *= 0.2;
    ::glLineWidth(3);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1,1,1);
    ::glBegin(GL_LINES);
    ::glVertex3d(px,py,pz);
    ::glVertex3d(px+ex.x, py+ex.y, pz+ex.z);
    ::glEnd();
  }
  
  ::glEnable(GL_LIGHTING);
  {
    //    float gray[4] = {0.3,0.3,0.3,1};
    float gray[4] = {0.9f,0.9f,0.9f,1.f};    
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    float shine[4] = {0,0,0,0};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
    ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 127.0);
    //    ::glColor3d(1,1,1);
  }
  
  //  DrawTri3D_FaceNorm(aXYZ, aTri, aNorm);
  
  if( m_texName != 0 ){
    //    ::glColor3d(1.0, 0.8, 0.5);
    glEnable(GL_TEXTURE_2D);
    //	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    //	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE );
    //    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL );
    
    //glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
    //glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glShadeModel(GL_SMOOTH);
    glBindTexture(GL_TEXTURE_2D , m_texName);
  }
  
  //  ::glDisable(GL_TEXTURE_2D);
  
  //  ::glDisable(GL_LIGHTING);
  {
    float gray[4] = {1,1,0.5,1};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    float shine[4] = {1,1,1,1};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
    ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128.0);
    //    ::glColor3d(1,1,1);
  }
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itri=0;itri<aTri.size()/3;itri++){
    int i0 = aTri[itri*3+0];
    int i1 = aTri[itri*3+1];
    int i2 = aTri[itri*3+2];
//    const double* pn0 = aLocCoord.data()+i0*6;
//    const double* pn1 = aLocCoord.data()+i1*6;
//    const double* pn2 = aLocCoord.data()+i2*6;
    double r = 1.0;
    //
    ::glNormal3dv(aLocCoord.data()+i0*6);
    ::glTexCoord2d(aTex[i0*2+0]*r+0.5,aTex[i0*2+1]*r+0.5);
    //    ::glColor3d(pn0[0]*0.5+0.5, pn0[1]*0.5+0.5, pn0[2]*0.5+0.5);
    ::glVertex3dv(aXYZ.data()+i0*3);
    //
    ::glNormal3dv(aLocCoord.data()+i1*6);
    ::glTexCoord2d(aTex[i1*2+0]*r+0.5,aTex[i1*2+1]*r+0.5);
    //    ::glColor3d(pn1[0]*0.5+0.5, pn1[1]*0.5+0.5, pn1[2]*0.5+0.5);
    ::glVertex3dv(aXYZ.data()+i1*3);
    //
    ::glNormal3dv(aLocCoord.data()+i2*6);
    ::glTexCoord2d(aTex[i2*2+0]*r+0.5,aTex[i2*2+1]*r+0.5);
    //    ::glColor3d(pn2[0]*0.5+0.5, pn2[1]*0.5+0.5, pn2[2]*0.5+0.5);
    ::glVertex3dv(aXYZ.data()+i2*3);
  }
  ::glEnd();
  if( !is_lighting ){ ::glDisable(GL_LIGHTING); }
  ::glDisable(GL_TEXTURE_2D);
  
  if( is_lighting ){ ::glEnable( GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }
  if(  is_texture  ){ ::glEnable(GL_TEXTURE_2D); } 
}

int main(int argc,char* argv[])
{
  SetNewProblem();
  // ------------------------------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();

  m_texName = ReadPPM_SetTexture(std::string(PATH_INPUT_DIR)+"/dep.ppm");

  int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    if( iframe % 100 == 0 ){
      iker = (int)((aXYZ.size()/3.0)*(double)rand()/(1.0+RAND_MAX));
      MakeExpMap_Point
          (aTex,aLocCoord,
           iker,
           aXYZ,psup_ind,psup);
      iframe = 0;
    }
    iframe++;
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
