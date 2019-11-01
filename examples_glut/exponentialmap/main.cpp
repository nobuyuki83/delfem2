/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <queue>
#include "delfem2/msh.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vec3.h"

// --------

#if defined(__APPLE__) && (__MACH__)
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_color.h"
#include "delfem2/opengl/gl_tex.h"
#include "../glut_cam.h"

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
  CIndNode(int ino, double mesh_dist) :ino(ino), mesh_dist(mesh_dist){}
  bool operator < (const CIndNode& lhs) const {
    return this->mesh_dist < lhs.mesh_dist;
  }
public:
  int ino;
  double mesh_dist;
};

void GetLocalExpMap
(double& u, double& v,
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

void GetGeodesic
(double& u_pq, double& v_pq,
 // const double rot[4],
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

void MakeLocalCoord
(double lrx[3], double lry[3],
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

void MakeExpMap_Point
(std::vector<double>& aTex,
 std::vector<double>& aLocCoord,
 int iker,
 const std::vector<double>& aXYZ,
 std::vector<int>& psup_ind,
 std::vector<int>& psup)
{
  const unsigned int nXYZ = aXYZ.size()/3;
  std::vector<CNode> aNode;
  aNode.resize(nXYZ);
  for(unsigned int i=0;i<nXYZ*2;i++){ aTex[i]=0; }
  ////
  double LocCoord_ker[9];
  for(unsigned int i=0;i<6;i++){ LocCoord_ker[i] = aLocCoord[iker*6+i]; }
  Cross3D(LocCoord_ker+6, LocCoord_ker+0, LocCoord_ker+3);
  ////
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
    ////
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
    ////
    for(int ipsup=psup_ind[ifix];ipsup<psup_ind[ifix+1];++ipsup){
      unsigned int ino1 = psup[ipsup];
      if( aNode[ino1].itype == 2 ) continue;
      const double len = Distance3D(aXYZ.data()+ifix*3,aXYZ.data()+ino1*3);
      const double weight = 1.0/len;
      const double dist0 = aNode[ifix].mesh_dist+len;
      double u, v;
      GetGeodesic
      (u,v,
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
std::vector<int> psup_ind;
std::vector<int> psup;
int iker = 0;
CNav3D_GLUT nav;
bool is_animation = true;
bool is_lighting = false;
int m_texName;
    
// ---------------------------------------------------

void SetNewProblem()
{
  const unsigned int nprob = 2;
  static unsigned int iprob = 0;
  
  {
    Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_2k.ply", aXYZ, aTri);
    {
      double cx,cy,cz, wx,wy,wz;
      GetCenterWidth(cx,cy,cz,
                     wx,wy,wz,
                     aXYZ);
      Translate(-cx,-cy,-cz, aXYZ);
      double wm = wx;
      wm = ( wx > wm ) ? wx : wm;
      wm = ( wy > wm ) ? wy : wm;
      wm = ( wz > wm ) ? wz : wm;
      Scale(2.0/wm,aXYZ);
    }
  }
  {
    std::vector<int> elsup_ind, elsup;
    JArrayElemSurPoint_MeshElem(elsup_ind, elsup,
                             aTri.data(), aTri.size()/3, 3, aXYZ.size()/3);
    JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                                aTri.data(), elsup_ind, elsup, 3, aXYZ.size()/3);
  }
  {
    std::vector<double> aNorm(aXYZ.size());
    Normal_MeshTri3D(aNorm.data(),
                     aXYZ.data(), aXYZ.size()/3, aTri.data(),aTri.size()/3);
    const int np = aXYZ.size()/3;
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

void myGlutIdle(){
  glutPostRedisplay();
}

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  ::glMatrixMode(GL_PROJECTION);
  glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  nav.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y){
  nav.glutMouse(button, state, x, y);
}

void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key, x, y);
}

void myGlutDisplay(void)
{
  ::glClearColor(0.2f, .7f, .7f,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  DrawBackground();
  nav.SetGL_Camera();
  
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
    ::glTranslated(+px, +py, +pz);
    ::glutSolidSphere(0.05, 16, 16);
    ::glTranslated(-px, -py, -pz);
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
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE );
    //    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL );
    
    //glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
    //glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
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
    const double* pn0 = aLocCoord.data()+i0*6;
    const double* pn1 = aLocCoord.data()+i1*6;
    const double* pn2 = aLocCoord.data()+i2*6;
    double r = 1.0;
    ////
    ::glNormal3dv(aLocCoord.data()+i0*6);
    ::glTexCoord2d(aTex[i0*2+0]*r+0.5,aTex[i0*2+1]*r+0.5);
    //    ::glColor3d(pn0[0]*0.5+0.5, pn0[1]*0.5+0.5, pn0[2]*0.5+0.5);
    ::glVertex3dv(aXYZ.data()+i0*3);
    ////
    ::glNormal3dv(aLocCoord.data()+i1*6);
    ::glTexCoord2d(aTex[i1*2+0]*r+0.5,aTex[i1*2+1]*r+0.5);
    //    ::glColor3d(pn1[0]*0.5+0.5, pn1[1]*0.5+0.5, pn1[2]*0.5+0.5);
    ::glVertex3dv(aXYZ.data()+i1*3);
    ////
    ::glNormal3dv(aLocCoord.data()+i2*6);
    ::glTexCoord2d(aTex[i2*2+0]*r+0.5,aTex[i2*2+1]*r+0.5);
    //    ::glColor3d(pn2[0]*0.5+0.5, pn2[1]*0.5+0.5, pn2[2]*0.5+0.5);
    ::glVertex3dv(aXYZ.data()+i2*3);
  }
  ::glEnd();
  if( !is_lighting ){ ::glDisable(GL_LIGHTING); }
  ::glDisable(GL_TEXTURE_2D);
  
  /*
   ::glBegin(GL_TRIANGLES);
   for(int itri=0;itri<aETri.size();itri++){
   const int i0 = aETri[itri].v[0];
   const int i1 = aETri[itri].v[1];
   const int i2 = aETri[itri].v[2];
   MyGlNormal3dv(aEPo[i0].n); MyGlVertex3dv(aEPo[i0].p);
   MyGlNormal3dv(aEPo[i1].n); MyGlVertex3dv(aEPo[i1].p);
   MyGlNormal3dv(aEPo[i2].n); MyGlVertex3dv(aEPo[i2].p);
   }
   ::glEnd();        
   */
  /*
   ::glDisable(GL_LIGHTING);
   ::glColor3d(0,0,0);
   ::glBegin(GL_LINES);
   for(unsigned int itri=0;itri<aTri.size();itri++){  
   const unsigned int i1 = aTri[itri].v[0];
   const unsigned int i2 = aTri[itri].v[1];
   const unsigned int i3 = aTri[itri].v[2];
   MyGlVertex3dv(aPo[i1].p);     MyGlVertex3dv(aPo[i2].p); 
   MyGlVertex3dv(aPo[i2].p);     MyGlVertex3dv(aPo[i3].p); 
   MyGlVertex3dv(aPo[i3].p);     MyGlVertex3dv(aPo[i1].p); 
   }
   ::glEnd();        
   
   */
  
  
  if( is_lighting ){ ::glEnable( GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }
  if(  is_texture  ){ ::glEnable(GL_TEXTURE_2D); } 
  
  //	ShowFPS();
  /*
   {
   int viewport[4];
   ::glGetIntegerv(GL_VIEWPORT, viewport);
   int width = viewport[2];
   int height = viewport[3];
   ::glViewport(0,0, 300,300);
   ::glMatrixMode(GL_PROJECTION);
   ::glLoadIdentity();
   ::glMatrixMode(GL_MODELVIEW);
   ::glLoadIdentity();
   
   //    ::glDisable(GL_TEXTURE_2D);
   
   //    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
   //    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
   ::glDisable(GL_LIGHTING);
   ::glDisable(GL_TEXTURE_2D);
   
   //  ::glDisable(GL_LIGHTING);
   
   //    ::glColor3d(1,1,1);
   //    ::glBegin(GL_TRIANGLES);
   //    ::glVertex2d(-1,-1);
   //    ::glVertex2d(+1,-1);
   //    ::glVertex2d( 0,+1);
   //    ::glEnd();
   
   ::glBegin(GL_TRIANGLES);
   double R[9] = { 0,0,0, 0,0,0, 0,0,0 };
   for(int i=0;i<6;i++){ R[i] = aLocCoord[iker*6+i]; }
   Cross3D(R+6, R+0, R+3);
   {
   double q0[3];
   //      VecMat3(aLocCoord.data()+iker*6, R, q0);
   MatVec3D(R, aLocCoord.data()+iker*6, q0);
   }
   for(int itri=0;itri<aTri.size()/3;itri++){
   int i0 = aTri[itri*3+0];
   int i1 = aTri[itri*3+1];
   int i2 = aTri[itri*3+2];
   const double* pn0 = aLocCoord.data()+i0*6;
   const double* pn1 = aLocCoord.data()+i1*6;
   const double* pn2 = aLocCoord.data()+i2*6;
   //      const double* pn0 = aNormal.data()+i0*3;
   //      const double* pn1 = aNormal.data()+i1*3;
   //      const double* pn2 = aNormal.data()+i2*3;
   //      double pN0[3] = {pn0[0],pn0[1],pn0[2]};
   //      double pN1[3] = {pn1[0],pn1[1],pn1[2]};
   //      double pN2[3] = {pn2[0],pn2[1],pn2[2]};
   //      double pN0[3]; VecMat3(pn0, R, pN0);
   //      double pN1[3]; VecMat3(pn1, R, pN1);
   //      double pN2[3]; VecMat3(pn2, R, pN2);
   double pN0[3]; MatVec3D(R, pn0, pN0);
   double pN1[3]; MatVec3D(R, pn1, pN1);
   double pN2[3]; MatVec3D(R, pn2, pN2);
   const double* pp0 = aXYZ.data()+i0*3;
   const double* pp1 = aXYZ.data()+i1*3;
   const double* pp2 = aXYZ.data()+i2*3;
   double r = 1.0;
   double s = 2.0;
   if( aTex[i0*2+0] < -0.5 || aTex[i0*2+0] > 0.5 ) continue;
   if( aTex[i0*2+1] < -0.5 || aTex[i0*2+1] > 0.5 ) continue;
   if( aTex[i1*2+0] < -0.5 || aTex[i1*2+0] > 0.5 ) continue;
   if( aTex[i1*2+1] < -0.5 || aTex[i1*2+1] > 0.5 ) continue;
   if( aTex[i2*2+0] < -0.5 || aTex[i2*2+0] > 0.5 ) continue;
   if( aTex[i2*2+1] < -0.5 || aTex[i2*2+1] > 0.5 ) continue;
   //      ::glColor3d(pn0[0]*0.5+0.5, pn0[1]*0.5+0.5, pn0[2]*0.5+0.5);
   //      ::glTexCoord2d(aTex[i0*2+0]*r+0.5,aTex[i0*2+1]*r+0.5);
   //      ::glColor3dv(pn1);
   //      ::glColor3d(pn1[0]*0.5+0.5, pn1[1]*0.5+0.5, pn1[2]*0.5+0.5);
   //      ::glTexCoord2d(aTex[i1*2+0]*r+0.5,aTex[i1*2+1]*r+0.5);
   //      ::glColor3dv(pn2);
   //      ::glColor3d(pn2[0]*0.5+0.5, pn2[1]*0.5+0.5, pn2[2]*0.5+0.5);
   //      ::glTexCoord2d(aTex[i2*2+0]*r+0.5,aTex[i2*2+1]*r+0.5);
   ////
   ::glColor3d(pN0[1]*0.5+0.5, pN0[2]*0.5+0.5, pN0[0]*0.5+0.5);
   ::glVertex2d(aTex[i0*2+0]*s,aTex[i0*2+1]*s);
   ////
   ::glColor3d(pN1[1]*0.5+0.5, pN1[2]*0.5+0.5, pN1[0]*0.5+0.5);
   ::glVertex2d(aTex[i1*2+0]*s,aTex[i1*2+1]*s);
   ////
   ::glColor3d(pN2[1]*0.5+0.5, pN2[2]*0.5+0.5, pN2[0]*0.5+0.5);
   ::glVertex2d(aTex[i2*2+0]*s,aTex[i2*2+1]*s);
   }
   ::glEnd();
   
   
   ::glViewport(0,0,width,height);
   }
   */
  
  glutSwapBuffers();
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key){
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'a':
      if( is_animation ){ is_animation = false; }
      else{ is_animation = true; }
      break;
    case ' ':
    {
      iker = (int)((aXYZ.size()/3.0)*(double)rand()/(1.0+RAND_MAX));
      MakeExpMap_Point
      (aTex,aLocCoord,
       iker,
       aXYZ,psup_ind,psup);
    }
      break;
    case 'l':
      is_lighting = !is_lighting;
      if( is_lighting ){ 
        ::glEnable(GL_LIGHTING);
        std::cout << "put on light " << std::endl;
      }
      else{              
        ::glDisable(GL_LIGHTING);
      }
      break;
  }
  
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  // initialize glut
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("Discrete Exponential Map");
  
  // define call back functions
  glutIdleFunc(myGlutIdle);
  glutKeyboardFunc(myGlutKeyboard);
  glutDisplayFunc(myGlutDisplay);
  glutReshapeFunc(myGlutResize);
  glutSpecialFunc(myGlutSpecial);;
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  
  setSomeLighting();
  SetNewProblem();
  
  //  m_texName = ReadPPM_SetTexture("checkerboard.pnm");
  
  m_texName = ReadPPM_SetTexture(std::string(PATH_INPUT_DIR)+"/dep.ppm");
  
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  nav.camera.view_height = 1.5;
  
  glutMainLoop();
  return 0;
}
