
#pragma warning( disable : 4786 )
#define for if(0); else for

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <fstream>
#include <time.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "isosurface_stuffing.h"
#include "utility_glut.h"
#include "utility_glcolor.h"
#include "signed_distance_function.h"
#include "mesh_tri.h"

///////////////////////////////////////////////////////////////////////////////////////////


static bool IsAbovePlane(const double p[3], const double org[3], const double n[3])
{
  const double dot = (p[0]-org[0])*n[0] + (p[1]-org[1])*n[1] + (p[2]-org[2])*n[2];
  return dot > 0;
}

static inline void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const double invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

static void makeTetSurroundingRelationship
(std::vector<int>& aElSurRel,
 const std::vector<int>& aEl,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup)
{
  const int nfael = 4;
  const int npofa = 3;
  const int npoel = 4;
  const unsigned int noelFace[nfael][npofa] = {
    { 1, 2, 3 },
    { 0, 3, 2 },
    { 0, 1, 3 },
    { 0, 2, 1 } };
  
  const int np = (int)elsup_ind.size()-1;
  const int nEl = (int)aEl.size()/npoel;
  aElSurRel.assign(nEl*nfael*2,-1);
  
  std::vector<int> tmp_poin(np,0);
  unsigned int inpofa[npofa];
  for (int iel = 0; iel<nEl; iel++){
    for (int ifael = 0; ifael<nfael; ifael++){
      for (int ipofa = 0; ipofa<npofa; ipofa++){
        int int0 = noelFace[ifael][ipofa];
        inpofa[ipofa] = aEl[iel*npoel+int0];
        tmp_poin[inpofa[ipofa]] = 1;
      }
      const int ipoin0 = inpofa[0];
      bool iflg = false;
      for (int ielsup = elsup_ind[ipoin0]; ielsup<elsup_ind[ipoin0+1]; ielsup++){
        const int jelem0 = elsup[ielsup];
        if (jelem0==iel) continue;
        for (int jfael = 0; jfael<nfael; jfael++){
          iflg = true;
          for (int jpofa = 0; jpofa<npofa; jpofa++){
            int jnt0 = noelFace[jfael][jpofa];
            const int jpoin0 = aEl[jelem0*npoel+jnt0];
            if (tmp_poin[jpoin0]==0){ iflg = false; break; }
          }
          if (iflg){
            aElSurRel[iel*nfael*2+ifael*2+0] = jelem0;
            aElSurRel[iel*nfael*2+ifael*2+1] = jfael;
            break;
          }
        }
        if (iflg) break;
      }
      if (!iflg){
        aElSurRel[iel*nfael*2+ifael*2+0] = -1;
        aElSurRel[iel*nfael*2+ifael*2+1] = -1;
      }
      for (int ipofa = 0; ipofa<npofa; ipofa++){
        tmp_poin[inpofa[ipofa]] = 0;
      }
    }
  }
}
/*
static void makeElemSurroundingPoint
(std::vector<int>& elsup_ind,
 std::vector<int>& elsup,
 ////
 const std::vector<int>& aElem,
 int nPoEl,
 int nPo)
{
  const int nElem = (int)aElem.size()/nPoEl;
  elsup_ind.assign(nPo+1,0);
  for(int ielem=0;ielem<nElem;ielem++){
    for(int inoel=0;inoel<nPoEl;inoel++){
      int ino1 = aElem[ielem*nPoEl+inoel];
      if( ino1 == -1 ){ break; }
      elsup_ind[ino1+1] += 1;
    }
  }
  for(int ino=0;ino<nPo;++ino){
    elsup_ind[ino+1] += elsup_ind[ino];
  }
  int nelsup = elsup_ind[nPo];
  elsup.resize(nelsup);
  for(int ielem=0;ielem<nElem;ielem++){
    for(int inoel=0;inoel<nPoEl;inoel++){
      int ino1 = aElem[ielem*nPoEl+inoel];
      if( ino1 == -1 ){ break; }
      int ind1 = elsup_ind[ino1];
      elsup[ind1] = ielem;
      elsup_ind[ino1] += 1;
    }
  }
  for(int ino=nPo;ino>=1;ino--){
    elsup_ind[ino] = elsup_ind[ino-1];
  }
  elsup_ind[0] = 0;
}

static void makeTetSurroundingRelationship
(std::vector<int>& aTetSurRel,
 const std::vector<int>& aTet,
 const int nXYZ)
{
  std::vector<int> elsup_ind, elsup;
  makeElemSurroundingPoint(elsup_ind, elsup, aTet, 4, nXYZ);
  makeTetSurroundingRelationship(aTetSurRel, aTet, elsup_ind,elsup);
}
 */

void IsoSurfaceStuffing_DrawCut
(const std::vector<double>& aXYZ, const std::vector<int>& aTet, const std::vector<int>& aIsSurfNode,
 const double org[3], const double n[3])
{
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_DIFFUSE);
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itet=0;itet<(int)aTet.size()/4;itet++){
    const int ino0 = aTet[itet*4+0];
    const int ino1 = aTet[itet*4+1];
    const int ino2 = aTet[itet*4+2];
    const int ino3 = aTet[itet*4+3];
    const double p0[3] = {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]};
    const double p3[3] = {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]};
    if( IsAbovePlane(p0, org, n) ) continue;
    if( IsAbovePlane(p1, org, n) ) continue;
    if( IsAbovePlane(p2, org, n) ) continue;
    if( IsAbovePlane(p3, org, n) ) continue;
    double n[3], area;
    ////
    UnitNormalAreaTri3D(n, area, p0, p2, p1);
    if( aIsSurfNode[ino0]==1 && aIsSurfNode[ino2]==1 && aIsSurfNode[ino1]==1 ){ ::glColor3d(1,1,0); }
    else{ ::glColor3d(1,0,1); }
    ::glNormal3dv(n);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    /////
    UnitNormalAreaTri3D(n, area, p0, p1, p3);
    if( aIsSurfNode[ino0]==1 && aIsSurfNode[ino1]==1 && aIsSurfNode[ino3]==1 ){ ::glColor3d(1,1,0); }
    else{ ::glColor3d(1,0,1); }
    ::glNormal3dv(n);
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    /////
    UnitNormalAreaTri3D(n, area, p1, p2, p3);
    if( aIsSurfNode[ino1]==1 && aIsSurfNode[ino2]==1 && aIsSurfNode[ino3]==1 ){ ::glColor3d(1,1,0); }
    else{ ::glColor3d(1,0,1); }
    ::glNormal3dv(n);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
    /////
    UnitNormalAreaTri3D(n, area, p2, p0, p3);
    if( aIsSurfNode[ino2]==1 && aIsSurfNode[ino0]==1 && aIsSurfNode[ino3]==1 ){ ::glColor3d(1,1,0); }
    else{ ::glColor3d(1,0,1); }
    ::glNormal3dv(n);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  bool is_lighting = glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itet=0;itet<aTet.size()/4;itet++){
    const int ino0 = aTet[itet*4+0];
    const int ino1 = aTet[itet*4+1];
    const int ino2 = aTet[itet*4+2];
    const int ino3 = aTet[itet*4+3];
    const double p0[3] = {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]};
    const double p3[3] = {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]};
    if( IsAbovePlane(p0, org, n) ) continue;
    if( IsAbovePlane(p1, org, n) ) continue;
    if( IsAbovePlane(p2, org, n) ) continue;
    if( IsAbovePlane(p3, org, n) ) continue;
    ////
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  ////
  /*
   ::glColor3d(0,0,0);
   ::glPointSize(3);
   ::glBegin(GL_POINTS);
   for(unsigned int ino=0;ino<nXYZ_;ino++){
   ::glVertex3dv(pXYZ_+ino*3);
   }
   ::glEnd();
   */
  if( is_lighting ){ glEnable(GL_LIGHTING); }
}


void IsoSurfaceStuffing_DrawCut2
(const std::vector<double>& aXYZ, const std::vector<int>& aTet, const std::vector<CColor>& aColor,
 const double org[3], const double n[3])
{
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_DIFFUSE);
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itet=0;itet<(int)aTet.size()/4;itet++){
    const int ino0 = aTet[itet*4+0];
    const int ino1 = aTet[itet*4+1];
    const int ino2 = aTet[itet*4+2];
    const int ino3 = aTet[itet*4+3];
    const double p0[3] = {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]};
    const double p3[3] = {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]};
    if( IsAbovePlane(p0, org, n) ) continue;
    if( IsAbovePlane(p1, org, n) ) continue;
    if( IsAbovePlane(p2, org, n) ) continue;
    if( IsAbovePlane(p3, org, n) ) continue;
//    ::glColor3d(1,1,0);
    aColor[itet].glColorDiffuse();
    ////
    double n[3], area;
    UnitNormalAreaTri3D(n, area, p0, p2, p1);
    ::glNormal3dv(n);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    /////
    UnitNormalAreaTri3D(n, area, p0, p1, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    /////
    UnitNormalAreaTri3D(n, area, p1, p2, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
    /////
    UnitNormalAreaTri3D(n, area, p2, p0, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  bool is_lighting = glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itet=0;itet<aTet.size()/4;itet++){
    const int ino0 = aTet[itet*4+0];
    const int ino1 = aTet[itet*4+1];
    const int ino2 = aTet[itet*4+2];
    const int ino3 = aTet[itet*4+3];
    const double p0[3] = {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]};
    const double p3[3] = {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]};
    if( IsAbovePlane(p0, org, n) ) continue;
    if( IsAbovePlane(p1, org, n) ) continue;
    if( IsAbovePlane(p2, org, n) ) continue;
    if( IsAbovePlane(p3, org, n) ) continue;
    ////
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  ////
  /*
   ::glColor3d(0,0,0);
   ::glPointSize(3);
   ::glBegin(GL_POINTS);
   for(unsigned int ino=0;ino<nXYZ_;ino++){
   ::glVertex3dv(pXYZ_+ino*3);
   }
   ::glEnd();
   */
  if( is_lighting ){ glEnable(GL_LIGHTING); }
}

CGlutWindowManager win;
double vis_cut_org[3] = {-0.0, 0.0, 0.0};
double vis_cut_nrm[3] = {-0.8, 0.1, 0.0};
//double vis_cut_time = 0;
bool is_animation = false;
double cur_time = 0;

int imode_draw = 1;

std::vector<double> aXYZ;
std::vector<int> aTet;
std::vector<int> aTetSurRel;
std::vector<int> aIsOnSurfXYZ;
std::vector<CColor> aColor;

static inline double TetVolume3D
(const double v1[3], const double v2[3], const double v3[3], const double v4[3] )
{
  return
  (   ( v2[0] - v1[0] )*( ( v3[1] - v1[1] )*( v4[2] - v1[2] ) - ( v4[1] - v1[1] )*( v3[2] - v1[2] ) )
   -	( v2[1] - v1[1] )*( ( v3[0] - v1[0] )*( v4[2] - v1[2] ) - ( v4[0] - v1[0] )*( v3[2] - v1[2] ) )
   +	( v2[2] - v1[2] )*( ( v3[0] - v1[0] )*( v4[1] - v1[1] ) - ( v4[0] - v1[0] )*( v3[1] - v1[1] ) )
   ) * 0.16666666666666666666666666666667;
}

void SetProblem()
{
  const unsigned int nprob = 6;
  static int iprob = 0;
  
  if( iprob == 0 ){
    class CSphere : public CInputIsosurfaceStuffing
    {
    public:
      CSphere(double rad){
        sp.radius_ = rad;
      }
      virtual double SignedDistance(double px, double py, double pz) const{
        double n[3];
        return sp.Projection(px,py,pz,n);
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                        double px, double py, double pz) const
      {
        sdf = this->SignedDistance(px,py,pz);
        const double rad0 = sp.radius_;
        const double rad1 = sqrt(px*px+py*py+pz*pz);
        if( rad1 > rad0*0.8 ){ ilevel_vol = 0; }
        else{ ilevel_vol = 2; }
        ilevel_srf = -1;
        nlayer = 1;
      }
    public:
      CSignedDistanceField3D_Sphere sp;
    } input(0.8);
    double elen = 0.25;
    int ndiv = 10;
    double cent[3] = {0,0,0};
    double org[3] = {cent[0]-ndiv*elen*0.5, cent[1]-ndiv*elen*0.5, cent[2]-ndiv*elen*0.5};
    std::vector<int> aTetLattice;
    std::vector<CPointLattice> aPointLattice;
    makeBackgroundLattice(aPointLattice,aTetLattice,
                          input, elen, ndiv, org);
    aTet = aTetLattice;
    aXYZ.resize(aPointLattice.size()*3);
    for(int ip=0;ip<aPointLattice.size();++ip){
      aXYZ[ip*3+0] = aPointLattice[ip].pos[0];
      aXYZ[ip*3+1] = aPointLattice[ip].pos[1];
      aXYZ[ip*3+2] = aPointLattice[ip].pos[2];
    }
    
    double tot_vol = 0;
    for(int itet=0;itet<aTet.size()/4;itet++){
      const int i0 = aTet[itet*4+0];
      const int i1 = aTet[itet*4+1];
      const int i2 = aTet[itet*4+2];
      const int i3 = aTet[itet*4+3];
      const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
      const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
      const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
      const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
      double vol = TetVolume3D(p0, p1, p2, p3);
      if( vol < 0 ){
        std::cout << "nega : " << itet << " " << vol << std::endl;
      }
      tot_vol += vol;
    }
    
    std::cout << " vaolume comp " << (ndiv*elen)*(ndiv*elen)*(ndiv*elen) << " " << tot_vol << std::endl;
    std::cout << "unmatched face should be: " << (ndiv*ndiv*2)*6 << std::endl;
  }
  else if( iprob == 1 ){
    class CSphere : public CInputIsosurfaceStuffing
    {
    public:
      CSphere(double rad){
        sp.radius_ = rad;
      }
      virtual double SignedDistance(double px, double py, double pz) const{
        double n[3];
        return sp.Projection(px,py,pz,n);
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                         double px, double py, double pz) const
      {
        sdf = SignedDistance(px, py, pz);
        ilevel_vol = -1;
        ilevel_srf = 2;
        nlayer = 2;
      }
    public:
      CSignedDistanceField3D_Sphere sp;
    };
    double rad = 1.5;
    CSphere sphere(rad);
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                       sphere, 0.5, rad*4.0, cent);
  }
  else if( iprob == 2 ){
    class CSphere : public CInputIsosurfaceStuffing
    {
    public:
      CSphere(double rad){
        sp.radius_ = rad;
      }
      virtual double SignedDistance(double px, double py, double pz) const{
        double n[3];
        return sp.Projection(px,py,pz,n);
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                         double px, double py, double pz) const
      {
        sdf = SignedDistance(px, py, pz);
        if( sdf > sp.radius_*0.5 ){ ilevel_vol = 0; }
        else{ ilevel_vol = 2; }
        ilevel_srf = -1;
        nlayer = 1;
      }
    public:
      CSignedDistanceField3D_Sphere sp;
    };
    double rad = 1.5;
    CSphere sphere(rad);
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                       sphere, 0.5, rad*4.0, cent);
  }
  else if( iprob ==  3 ){
    class CBox : public CInputIsosurfaceStuffing
    {
    public:
      CBox(double hwx, double hwy, double hwz){
        bx.hwx = hwx;
        bx.hwy = hwy;
        bx.hwz = hwz;
      }
      virtual double SignedDistance(double px, double py, double pz) const{
        double n[3];
        return bx.Projection(px,py,pz,n);
      }
    public:
      CSignedDistanceField3D_Box bx;
    };
    const double hwx = 0.91;
    const double hwy = 0.61;
    const double hwz = 0.41;
    double cent[3] = {0,0,0};
    CBox box(hwx,hwy,hwz);
    IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                       box, 0.1, 2.00, cent);
  }
  else if( iprob == 4 ){
    class CCavSphere : public CInputIsosurfaceStuffing
    {
    public:
      CCavSphere(){
        const double hwx = 0.91;
        const double hwy = 0.61;
        const double hwz = 0.61;
        box.hwx = hwx;
        box.hwy = hwy;
        box.hwz = hwz;
        ////
        const double rad = 0.2;
        sphere.radius_ = rad;
      }
      virtual double SignedDistance(double x, double y, double z ) const {
        double n[3];
        double dist0 = -sphere.Projection(x, y, z, n);
        double cx = 0.5;
        double cy = 0.0;
        double cz = 0.0;
        double dist1 = box.Projection(x-cx, y-cy, z-cz, n);
        return (dist0<dist1) ? dist0 : dist1;
      }
    public:
      CSignedDistanceField3D_Box box;
      CSignedDistanceField3D_Sphere sphere;
    } cav_sphere;
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet, aIsOnSurfXYZ,
                       cav_sphere, 0.1, 3.00, cent);
  }
  if( iprob == 5 ){
    class CMesh : public CInputIsosurfaceStuffing
    {
    public:
      virtual double SignedDistance(double x, double y, double z) const {
        double n[3];
        return sdf_mesh.Projection(x, y, z,n);
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                         double px, double py, double pz) const
      {
        sdf = this->SignedDistance(px,py,pz);
        ilevel_vol = 0;
        ilevel_srf = 2;
        nlayer = 2;
      }
    public:
      CSignedDistanceField3D_Mesh sdf_mesh;
    } mesh;
    {
      std::vector<int> aTri;
      std::vector<double> aXYZ_Tri;
      Load_Ply("bunny_1k.ply", aXYZ_Tri, aTri);
      Normalize(aXYZ_Tri,2.3);
      mesh.sdf_mesh.SetMesh(aTri, aXYZ_Tri);
      mesh.sdf_mesh.BuildBoxel();
    }
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet, aIsOnSurfXYZ,
                       mesh, 0.18, 3.0, cent);
  }
  
  makeTetSurroundingRelationship(aTetSurRel, aTet, (int)aXYZ.size()/3);
  
  aColor.resize(aTet.size()/4);
  for(int it=0;it<aTet.size()/4;++it){
    aColor[it].setRandomVividColor();
  }
  
  //  for(int ixyz=0;ixyz<aXYZ.size()/3;++ixyz){
  //    std::cout << ixyz << " " << aXYZ[ixyz*3+0] << " " << aXYZ[ixyz*3+1] << " " << aXYZ[ixyz*3+2] << std::endl;
  //  }
  {
    int icnt = 0;
    for(int itet=0;itet<aTet.size()/4;++itet){
      for(int iface=0;iface<4;++iface){
        if( aTetSurRel[itet*8+iface*2+0]==-1){
          icnt++;
        }
      }
    }
    std::cout << "unmatch face" << icnt << std::endl;
  }
  
  
  iprob++;
  if( iprob == nprob ){ iprob = 0; }
}


void DrawTetMesh()
{
  ::glColor3d(0,0,0);
  ::glPointSize(3);
  ::glBegin(GL_POINTS);
  for(unsigned int ino=0;ino<aXYZ.size()/4;ino++){
    const double p[3] = {aXYZ[ino*3+0], aXYZ[ino*3+1], aXYZ[ino*3+2]};
    ::glVertex3dv(p);
  }
  ::glEnd();
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itet=0;itet<aTet.size()/4;itet++){
    const int ino0 = aTet[itet*4+0];
    const int ino1 = aTet[itet*4+1];
    const int ino2 = aTet[itet*4+2];
    const int ino3 = aTet[itet*4+3];
    const double p0[3] = {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]};
    const double p3[3] = {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]};
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    ////
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    ////
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
    ////
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itet=0;itet<aTet.size()/4;itet++){
    const int ino0 = aTet[itet*4+0];
    const int ino1 = aTet[itet*4+1];
    const int ino2 = aTet[itet*4+2];
    const int ino3 = aTet[itet*4+3];
    const double p0[3] = {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]};
    const double p3[3] = {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]};
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
  }
  ::glEnd();
}

/*
void DrawFrame()
{
  ::glDisable(GL_LIGHTING);
  ::glPointSize(3);
  ::glBegin(GL_POINTS);
  for(int ip=0;ip<(int)aPoint.size();++ip){
    ::glColor3d(0,0,0);
    ::glVertex3dv(aPoint[ip].pos);
  }
  ::glEnd();
  ////
  ::glPointSize(5);
  ::glBegin(GL_POINTS);
  for(int ic=0;ic<(int)aCell.size();++ic){
    const CCell& c = aCell[ic];
    int ip = c.aIP[26];
//    if( c.is_leaf == true ){ ::glColor3d(1,0,0); }
//    else{                    ::glColor3d(0,0,0); }
    ::glColor3d(0,0,0);
    ::glVertex3dv(aPoint[ip].pos);
  }
  ::glEnd();
  ////
  for(int ic=0;ic<(int)aCell.size();++ic){
    const CCell& c = aCell[ic];
    //    if( c.ilevel == 0 ) continue;
    {
      int ip0 = c.aIP[0];
      int ip1 = c.aIP[1];
      int ip2 = c.aIP[2];
      int ip3 = c.aIP[3];
      int ip4 = c.aIP[4];
      int ip5 = c.aIP[5];
      int ip6 = c.aIP[6];
      int ip7 = c.aIP[7];
      
      ::glColor3d(0,0,0);
      ::glLineWidth(1);
      ::glBegin(GL_LINES);
      if( ip0>=0 && ip1>=0 ){ ::glVertex3dv(aPoint[ip0].pos); ::glVertex3dv(aPoint[ip1].pos); }
      if( ip1>=0 && ip3>=0 ){ ::glVertex3dv(aPoint[ip1].pos); ::glVertex3dv(aPoint[ip3].pos); }
      if( ip3>=0 && ip2>=0 ){ ::glVertex3dv(aPoint[ip3].pos); ::glVertex3dv(aPoint[ip2].pos); }
      if( ip2>=0 && ip0>=0 ){ ::glVertex3dv(aPoint[ip2].pos); ::glVertex3dv(aPoint[ip0].pos); }
      ////
      if( ip4>=0 && ip5>=0 ){ ::glVertex3dv(aPoint[ip4].pos); ::glVertex3dv(aPoint[ip5].pos); }
      if( ip5>=0 && ip7>=0 ){ ::glVertex3dv(aPoint[ip5].pos); ::glVertex3dv(aPoint[ip7].pos); }
      if( ip7>=0 && ip6>=0 ){ ::glVertex3dv(aPoint[ip7].pos); ::glVertex3dv(aPoint[ip6].pos); }
      if( ip6>=0 && ip4>=0 ){ ::glVertex3dv(aPoint[ip6].pos); ::glVertex3dv(aPoint[ip4].pos); }
      ////
      if( ip0>=0 && ip4>=0 ){ ::glVertex3dv(aPoint[ip0].pos); ::glVertex3dv(aPoint[ip4].pos); }
      if( ip1>=0 && ip5>=0 ){ ::glVertex3dv(aPoint[ip1].pos); ::glVertex3dv(aPoint[ip5].pos); }
      if( ip2>=0 && ip6>=0 ){ ::glVertex3dv(aPoint[ip2].pos); ::glVertex3dv(aPoint[ip6].pos); }
      if( ip3>=0 && ip7>=0 ){ ::glVertex3dv(aPoint[ip3].pos); ::glVertex3dv(aPoint[ip7].pos); }
      ::glEnd();
    }
    
    for(int iica=0;iica<6;++iica){
      int ica0 = c.aIC_Adj[iica];
      if( ica0 < 0 ){ continue; }
      const CCell& ca = aCell[ica0];
      if( ca.ilevel != c.ilevel ) continue;
//      if( !ca.is_leaf || !c.is_leaf ) continue;
      const int ip = c.aIP[26];
      int iq = ca.aIP[26];
      ::glLineWidth(1);
      ::glBegin(GL_LINES);
      ::glColor3d(0,1,0);
      ::glVertex3dv(aPoint[ip].pos);
      ::glVertex3dv(aPoint[iq].pos);
      ::glEnd();
    }
    
    
    {
      const int ip = c.aIP[26];
      for(int iiq=0;iiq<8;++iiq){
        int iq = c.aIP[iiq];
        if( iq == -1 ) continue;
        ::glLineWidth(1);
        ::glBegin(GL_LINES);
        ::glColor3d(1,0,0);
        ::glVertex3dv(aPoint[ip].pos);
        ::glVertex3dv(aPoint[iq].pos);
        ::glEnd();
      }
    }
  }
}
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void myGlutResize(int w, int h)
{
  ::glViewport(0, 0, w, h);
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  win.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y){
  win.glutMouse(button,state,x,y);
}

void SetProblem();
void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
		case 'q':
		case 'Q':
		case '\033':
			exit(0);  /* '\033' ? ESC ? ASCII ??? */
		case 'a':
			is_animation = !is_animation;
			break;
    case 'd':
      imode_draw = (imode_draw+1)%3;
      break;
		case 'b':
			cur_time = 0;
			break;
		case ' ':	
			SetProblem();
			break;
	}
}


void myGlutSpecial(int Key, int x, int y)
{
  win.glutSpecial(Key,x,y);
}

void myGlutIdle(){
  if( is_animation ){
    cur_time += 0.005;
    if( cur_time > 1 ){ cur_time = 0.0; }
    vis_cut_org[0] = -1*(1-cur_time) + 1*cur_time;
  }
	::glutPostRedisplay();
}

////////////////

void myGlutDisplay(void)
{
//	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);  
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);
	
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
  win.SetGL_Camera();


  if( imode_draw == 0 ){
    ::glDisable(GL_LIGHTING);
//    DrawFrame();
  }
  else if( imode_draw == 1 ){
    ::glEnable(GL_LIGHTING);
    IsoSurfaceStuffing_DrawCut2(aXYZ,aTet,aColor,
                               vis_cut_org, vis_cut_nrm);
  }
  else if( imode_draw == 2 ){
    const int noelTetFace[4][3] = {
      { 1, 2, 3 },
      { 0, 3, 2 },
      { 0, 1, 3 },
      { 0, 2, 1 },
    };
    ::glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    ::glBegin(GL_TRIANGLES);
    for(int itet=0;itet<(int)aTet.size()/4;++itet){
      const int ino0 = aTet[itet*4+0];
      const int ino1 = aTet[itet*4+1];
      const int ino2 = aTet[itet*4+2];
      const int ino3 = aTet[itet*4+3];
      const double aP[4][3] = {
        {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]},
        {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]},
        {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]},
        {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]} };
      if( IsAbovePlane(aP[0], vis_cut_org, vis_cut_nrm) ) continue;
      if( IsAbovePlane(aP[1], vis_cut_org, vis_cut_nrm) ) continue;
      if( IsAbovePlane(aP[2], vis_cut_org, vis_cut_nrm) ) continue;
      if( IsAbovePlane(aP[3], vis_cut_org, vis_cut_nrm) ) continue;
      aColor[itet].glColorDiffuse();
      for(int iface=0;iface<4;++iface){
        const double* p0 = aP[ noelTetFace[iface][0] ];
        const double* p1 = aP[ noelTetFace[iface][1] ];
        const double* p2 = aP[ noelTetFace[iface][2] ];
        if( aTetSurRel[itet*8+iface*2+0] == -1 ){
          double n[3], area;
          UnitNormalAreaTri3D(n, area, p0, p1, p2);
          ::glNormal3dv(n);
          ::glVertex3dv(p0);
          ::glVertex3dv(p2);
          ::glVertex3dv(p1);
        }
      }
    }
    ::glEnd();
    
  }
	
//	ShowFPS();
	::glutSwapBuffers();
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{	
	// Initialize GLUT
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");
	
	// Setting call back function
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);
	
	SetProblem();
  win.camera.view_height = 2;
  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  ////
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_LIGHT0);            
  float light0pos[4] = {0.5,0.5,+20,0};
  ::glLightfv(GL_LIGHT0, GL_POSITION, light0pos);      
  float white[3] = {1.0,1.0,1.0};
  ::glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
  
  ::glEnable(GL_LIGHT1);    
  float light1pos[4] = {0.5,0.5,-20,0}; 
  ::glLightfv(GL_LIGHT1, GL_POSITION, light1pos);          
  float white2[3] = {1.0,1.0,0.0};
  ::glLightfv(GL_LIGHT1, GL_DIFFUSE, white2);        
  ////
  
  
	glutMainLoop();
	return 0;
}
