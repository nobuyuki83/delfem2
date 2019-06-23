/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <map>

#include "delfem2/mshsrch_v3.h"

CVector3 CPointElemSolid::getPos_Tet
(const std::vector<double> &aXYZ,
 const std::vector<int> &aTet) const
{
  assert(ielem>=0&&ielem<(int)aTet.size()/4);
  int ip0 = aTet[ielem*4+0];
  int ip1 = aTet[ielem*4+1];
  int ip2 = aTet[ielem*4+2];
  int ip3 = aTet[ielem*4+3];
  const CVector3 p0(aXYZ[ip0*3+0], aXYZ[ip0*3+1], aXYZ[ip0*3+2]);
  const CVector3 p1(aXYZ[ip1*3+0], aXYZ[ip1*3+1], aXYZ[ip1*3+2]);
  const CVector3 p2(aXYZ[ip2*3+0], aXYZ[ip2*3+1], aXYZ[ip2*3+2]);
  const CVector3 p3(aXYZ[ip3*3+0], aXYZ[ip3*3+1], aXYZ[ip3*3+2]);
  return r0*p0+r1*p1+r2*p2+(1.0-r0-r1-r2)*p3;
}

void CPointElemSolid::setPos_Tet
(int it0,
 const CVector3 &q,
 const std::vector<double> &aXYZ,
 const std::vector<int> &aTet)
{
  assert(it0>=0&&it0<(int)aTet.size()/4);
  int ip0 = aTet[it0*4+0];
  int ip1 = aTet[it0*4+1];
  int ip2 = aTet[it0*4+2];
  int ip3 = aTet[it0*4+3];
  const CVector3 p0(aXYZ[ip0*3+0], aXYZ[ip0*3+1], aXYZ[ip0*3+2]);
  const CVector3 p1(aXYZ[ip1*3+0], aXYZ[ip1*3+1], aXYZ[ip1*3+2]);
  const CVector3 p2(aXYZ[ip2*3+0], aXYZ[ip2*3+1], aXYZ[ip2*3+2]);
  const CVector3 p3(aXYZ[ip3*3+0], aXYZ[ip3*3+1], aXYZ[ip3*3+2]);
  double v0 = volume_Tet( q,p1,p2,p3);
  double v1 = volume_Tet(p0, q,p2,p3);
  double v2 = volume_Tet(p0,p1, q,p3);
  //    double v3 = volume_Tet(p0,p1,p2, q);
  double vt = volume_Tet(p0,p1,p2,p3);
  this->ielem = it0;
  this->r0 = v0/vt;
  this->r1 = v1/vt;
  this->r2 = v2/vt;
}

CVector3 CPointElemSurf::getPos_Tri(const std::vector<double>& aXYZ, const std::vector<int>& aTri) const
{
  assert(itri>=0&&itri<(int)aTri.size()/3);
  int ip0 = aTri[itri*3+0];
  int ip1 = aTri[itri*3+1];
  int ip2 = aTri[itri*3+2];
  const CVector3 p0(aXYZ[ip0*3+0], aXYZ[ip0*3+1], aXYZ[ip0*3+2]);
  const CVector3 p1(aXYZ[ip1*3+0], aXYZ[ip1*3+1], aXYZ[ip1*3+2]);
  const CVector3 p2(aXYZ[ip2*3+0], aXYZ[ip2*3+1], aXYZ[ip2*3+2]);
  return r0*p0+r1*p1+(1.0-r0-r1)*p2;
}

CVector3 CPointElemSurf::getPos_TetFace
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTet,
 const std::vector<int>& aTetFace) const
{
  const int noelTetFace[4][3] ={
    { 1, 2, 3 },
    { 0, 3, 2 },
    { 0, 1, 3 },
    { 0, 2, 1 } };
  int itet = aTetFace[itri*2+0];
  int iface = aTetFace[itri*2+1];
  double r2 = (1-r0-r1);
  int ielno0 = noelTetFace[iface][0];
  int ielno1 = noelTetFace[iface][1];
  int ielno2 = noelTetFace[iface][2];
  int iq0 = aTet[itet*4+ielno0];
  int iq1 = aTet[itet*4+ielno1];
  int iq2 = aTet[itet*4+ielno2];
  CVector3 p;
  p.x = r0*aXYZ[iq0*3+0]+r1*aXYZ[iq1*3+0]+r2*aXYZ[iq2*3+0];
  p.y = r0*aXYZ[iq0*3+1]+r1*aXYZ[iq1*3+1]+r2*aXYZ[iq2*3+1];
  p.z = r0*aXYZ[iq0*3+2]+r1*aXYZ[iq1*3+2]+r2*aXYZ[iq2*3+2];
  return p;
}









/*
 void FindValueInTet(std::vector<CPointTet>& mapXYZ2Tet,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aXYZTet,
 const std::vector<int>& aTet)
 {
 const int noelTetFace[4][3] ={
 { 1, 2, 3 },
 { 0, 3, 2 },
 { 0, 1, 3 },
 { 0, 2, 1 } };
 
 std::vector<int> aTetSurRel;
 makeSurroundingRelationship(aTetSurRel,
 aTet,FEMELEM_TET,
 (int)aXYZTet.size()/3);
 
 std::vector<int> aTetFaceSrf;
 {
 aTetFaceSrf.clear();
 for(int itet=0;itet<aTet.size()/4;++itet){
 for(int iface=0;iface<4;++iface){
 if( aTetSurRel[itet*8+iface*2+0] != -1 ) continue;
 aTetFaceSrf.push_back(itet);
 aTetFaceSrf.push_back(iface);
 }
 }
 }
 
 {
 const int np = aXYZ.size()/3;
 mapXYZ2Tet.assign(np,CPointTet());
 for(int ip=0;ip<np;++ip){
 const CVector3 p0(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
 if( ip!=0 ){
 mapXYZ2Tet[ip] = findPointInTetMesh(p0, mapXYZ2Tet[ip-1].itet, aXYZTet, aTet, aTetSurRel);
 }
 if( mapXYZ2Tet[ip].itet == -1 ){
 mapXYZ2Tet[ip] = findPointInTetMesh(p0, aXYZTet, aTet);
 }
 if( mapXYZ2Tet[ip].itet == -1 ){
 CPointTetFace ptf = findNearestTetFace(p0, aXYZTet, aTet, aTetFaceSrf);
 int itet = aTetFaceSrf[ptf.itetface*2+0];
 assert( itet>=0 && itet<(int)aTet.size()/4 );
 int iface = aTetFaceSrf[ptf.itetface*2+1];
 double aR[4] = {0,0,0,0};
 double r0 = ptf.r0;
 double r1 = ptf.r1;
 double r2 = 1-r0-r1;
 aR[noelTetFace[iface][0]] = r0;
 aR[noelTetFace[iface][1]] = r1;
 aR[noelTetFace[iface][2]] = r2;
 mapXYZ2Tet[ip].itet = itet;
 mapXYZ2Tet[ip].r0 = aR[0];
 mapXYZ2Tet[ip].r1 = aR[1];
 mapXYZ2Tet[ip].r2 = aR[2];
 }
 assert( mapXYZ2Tet[ip].itet>=0 && mapXYZ2Tet[ip].itet<(int)aTet.size()/4 );
 }
 std::cout << "end smpl" << std::endl;
 }
 }
 
 void FindValueInTri(std::vector<CPointTri>& mapXYZ2Tri,
 ////
 const std::vector<double>& aXYZ,
 const std::vector<double>& aXYZTri,
 const std::vector<int>& aTri)
 {
 std::vector<int> aTriSurRel;
 makeSurroundingRelationship(aTriSurRel,
 aTri,FEMELEM_TRI,
 (int)aXYZTri.size()/3);
 
 {
 const int np = aXYZ.size()/3;
 mapXYZ2Tri.assign(np,CPointTri());
 for(int ip=0;ip<np;++ip){
 const CVector3 p0(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
 if( ip!=0 ){
 //        mapXYZ2Tri[ip] = findPointInTriMesh(p0, mapXYZ2Tri[ip-1].itri, aXYZTri, aTri, aTriSurRel);
 }
 if( mapXYZ2Tri[ip].itri == -1 ){
 mapXYZ2Tri[ip] = findPointInTriMesh(p0, aXYZTri, aTri);
 }
 assert( mapXYZ2Tri[ip].itri>=0 && mapXYZ2Tri[ip].itri<(int)aTri.size()/3 );
 }
 std::cout << "end smpl tri" << std::endl;
 }
 }
 */


std::ostream &operator<<(std::ostream &output, const CPointElemSurf& v)
{
  output.setf(std::ios::scientific);
  output<<v.itri<<" "<<v.r0<<" "<<v.r1;
  return output;
}

std::istream &operator>>(std::istream &input, CPointElemSurf& v)
{
  input>>v.itri>>v.r0>>v.r1;
  return input;
}

/*
CVector3 MidPoint
(int itri,
 const std::vector<int>& aTri,
 const std::vector<double>& aXYZ)
{
  CVector3 p;
  int i0 = aTri[itri*3+0];
  int i1 = aTri[itri*3+1];
  int i2 = aTri[itri*3+2];
  p.x = (aXYZ[i0*3+0]+aXYZ[i1*3+0]+aXYZ[i2*3+0])/3.0;
  p.y = (aXYZ[i0*3+1]+aXYZ[i1*3+1]+aXYZ[i2*3+1])/3.0;
  p.z = (aXYZ[i0*3+2]+aXYZ[i1*3+2]+aXYZ[i2*3+2])/3.0;
  return p;
}

void weightInTriangle
(double& r0, double& r1,
 const CVector3& p,
 int itri,
 const std::vector<int>& aTri,
 const std::vector<double>& aXYZ)
{
  int ip0 = aTri[itri*3+0];
  int ip1 = aTri[itri*3+1];
  int ip2 = aTri[itri*3+2];
  const CVector3 p0(aXYZ[ip0*3+0], aXYZ[ip0*3+1], aXYZ[ip0*3+2]);
  const CVector3 p1(aXYZ[ip1*3+0], aXYZ[ip1*3+1], aXYZ[ip1*3+2]);
  const CVector3 p2(aXYZ[ip2*3+0], aXYZ[ip2*3+1], aXYZ[ip2*3+2]);
  double a0 = TriArea(p, p1, p2);
  double a1 = TriArea(p, p2, p0);
  double a2 = TriArea(p, p0, p1);
  r0 = a0/(a0+a1+a2);
  r1 = a1/(a0+a1+a2);
}
 */

bool intersectRay_Tri3D
(double& r0, double& r1,
 const CVector3& org, const CVector3& dir,
 const CVector3& p0, const CVector3& p1, const CVector3& p2)
{
  const double v0 = volume_Tet(p1, p2, org, org+dir);
  const double v1 = volume_Tet(p2, p0, org, org+dir);
  const double v2 = volume_Tet(p0, p1, org, org+dir);
  const double vt = v0+v1+v2;
  r0 = v0/vt;
  r1 = v1/vt;
  const double r2 = v2/vt;
  if( r0>0&&r1>0&&r2>0 ){
    return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

CPointElemSurf intersect_Ray_Tri3D
(double& depth,
 const CVector3& org, const CVector3& dir,
 int itri,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ)
{
  const int ip0 = aTri[itri*3+0];  assert(ip0>=0&&ip0<(int)aXYZ.size()/3);
  const int ip1 = aTri[itri*3+1];  assert(ip1>=0&&ip1<(int)aXYZ.size()/3);
  const int ip2 = aTri[itri*3+2];  assert(ip2>=0&&ip2<(int)aXYZ.size()/3);
  const CVector3 p0(aXYZ[ip0*3+0], aXYZ[ip0*3+1], aXYZ[ip0*3+2]);
  const CVector3 p1(aXYZ[ip1*3+0], aXYZ[ip1*3+1], aXYZ[ip1*3+2]);
  const CVector3 p2(aXYZ[ip2*3+0], aXYZ[ip2*3+1], aXYZ[ip2*3+2]);
  double r0, r1;
  bool res = intersectRay_Tri3D(r0,r1,
                                org, dir, p0,p1,p2);
  if( !res ) return CPointElemSurf(-1,0,0);
  double r2 = 1-r0-r1;
  CVector3 q0 = p0*r0+p1*r1+p2*r2;
  depth = -(q0-org)*dir/dir.DLength();
  return CPointElemSurf(itri,r0,r1);
}

CPointElemSurf intersect_Ray_MeshTriFlag3D
(const CVector3& org, const CVector3& dir,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ,
 int iflag,
 const std::vector<int>& aFlag)
{
  assert( aTri.size()/3 == aFlag.size() );
  std::map<double, CPointElemSurf> pickMap;
  for (int itri = 0; itri<(int)aTri.size()/3; itri++){
    if( aFlag[itri] != iflag ) continue;
    double depth;
    CPointElemSurf res = intersect_Ray_Tri3D(depth, org,dir, itri, aTri,aXYZ);
    if( res.itri<0 ){ continue; }
    pickMap.insert(std::make_pair(depth, res));
  }
  if (pickMap.empty()) return CPointElemSurf();
  return pickMap.begin()->second;
  /*
  int itri_pick = pickMap.begin()->second;
  double depth_pick = pickMap.begin()->first;
  p = -depth_pick*dir+org;
  return itri_pick;
   */
}

CPointElemSurf intersect_Ray_MeshTri3D
(const CVector3& org, const CVector3& dir,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ)
{
  std::map<double, CPointElemSurf> pickMap;
  for (int itri = 0; itri<(int)aTri.size()/3; itri++){
    double depth;
    CPointElemSurf res = intersect_Ray_Tri3D(depth, org,dir, itri, aTri,aXYZ);
    if( res.itri<0 ){ continue; }
    pickMap.insert(std::make_pair(depth, res));
  }
  if (pickMap.empty()) return CPointElemSurf(-1,0,0);
  return pickMap.begin()->second;
}

CPointElemSurf intersect_Ray_MeshTri3D
(const CVector3& org, const CVector3& dir,
 int itri_start, // starting triangle
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTriSurRel)
{
  int itri1 = itri_start;
  if (itri1<0||itri1>=(int)aTri.size()/3){ return CPointElemSurf(); }
  for (int itr = 0; itr<50; itr++){
    if (itri1==-1) return CPointElemSurf();
    int ip0 = aTri[itri1*3+0];
    int ip1 = aTri[itri1*3+1];
    int ip2 = aTri[itri1*3+2];
    const CVector3 p0(aXYZ[ip0*3+0], aXYZ[ip0*3+1], aXYZ[ip0*3+2]);
    const CVector3 p1(aXYZ[ip1*3+0], aXYZ[ip1*3+1], aXYZ[ip1*3+2]);
    const CVector3 p2(aXYZ[ip2*3+0], aXYZ[ip2*3+1], aXYZ[ip2*3+2]);
    const double v0 = volume_Tet(p1, p2, org, org+dir);
    const double v1 = volume_Tet(p2, p0, org, org+dir);
    const double v2 = volume_Tet(p0, p1, org, org+dir);
    if (v0>0&&v1>0&&v2>0){
      double r0 = v0/(v0+v1+v2);
      double r1 = v1/(v0+v1+v2);
//      double r2 = v2/(v0+v1+v2);
      return CPointElemSurf(itri1,r0,r1);
    }
    if (v0<v1 && v0<v2){      itri1 = aTriSurRel[itri1*6+0*2+0]; }
    else if (v1<v0 && v1<v2){ itri1 = aTriSurRel[itri1*6+1*2+0]; }
    else{                     itri1 = aTriSurRel[itri1*6+2*2+0]; }
  }
  ////
  return CPointElemSurf();
}

CPointElemSurf intersect_Ray_MeshTri3D
(const CVector3& org, const CVector3& dir,
 const CPointElemSurf& ptri0,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTriSurRel)
{
  CPointElemSurf ptri;
  if( ptri0.itri != -1 ){
    ptri = intersect_Ray_MeshTri3D(org,dir, ptri0.itri,aTri,aXYZ,aTriSurRel);
  }
  if( ptri.itri == -1 ){
    ptri = intersect_Ray_MeshTri3D(org,dir,aTri,aXYZ);
  }
  if( ptri.itri == -1 ) return ptri;
  return ptri;
}


//////////////////////////////////////////////////////////////////////////////////////////

CPointElemSurf nearest_Point_MeshTri3D
(const CVector3& q,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri)
{
  CPointElemSurf pes;
  double min_dist = -1;
  for(int it=0;it<aTri.size()/3;++it){
    int ip0 = aTri[it*3+0];
    int ip1 = aTri[it*3+1];
    int ip2 = aTri[it*3+2];
    const CVector3 p0 = CVector3(aXYZ[ip0*3+0], aXYZ[ip0*3+1], aXYZ[ip0*3+2])-q;
    const CVector3 p1 = CVector3(aXYZ[ip1*3+0], aXYZ[ip1*3+1], aXYZ[ip1*3+2])-q;
    const CVector3 p2 = CVector3(aXYZ[ip2*3+0], aXYZ[ip2*3+1], aXYZ[ip2*3+2])-q;
    double r0,r1;
    CVector3 p_min = nearest_Origin_Tri(r0,r1, p0,p1,p2);
    double dist = (p_min-q).DLength();
    if( min_dist<0 || dist < min_dist ){
      min_dist = dist;
      pes = CPointElemSurf(it,r0,r1);
    }
  }
  assert( pes.itri != -1 );
  /////
  return pes;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

CPointElemSolid nearest_Point_MeshTet3D
(const CVector3& q,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTet)
{
  const double eps = 1.0e-4;
  const int ntet = (int)aTet.size()/4;
  for(int itet=0;itet<ntet;++itet){
    CPointElemSolid pt;
    pt.setPos_Tet(itet,q,aXYZ,aTet);
    if( pt.isInside(-eps) ){ return pt; }
  }
  return CPointElemSolid();
}

CPointElemSolid nearest_Point_MeshTet3D
(const CVector3& p,
 int itet_start, // starting triangle
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTet,
 const std::vector<int>& aTetSurRel)
{
  const double eps = 1.0e-4;
  int itet1 = itet_start;
  if (itet1<0||itet1>=(int)aTet.size()/4){ return CPointElemSolid(); }
  for (int itr = 0; itr<50; itr++){
    if (itet1==-1) return CPointElemSolid();
    int ip0 = aTet[itet1*4+0];
    int ip1 = aTet[itet1*4+1];
    int ip2 = aTet[itet1*4+2];
    int ip3 = aTet[itet1*4+3];
    const CVector3 p0(aXYZ[ip0*3+0], aXYZ[ip0*3+1], aXYZ[ip0*3+2]);
    const CVector3 p1(aXYZ[ip1*3+0], aXYZ[ip1*3+1], aXYZ[ip1*3+2]);
    const CVector3 p2(aXYZ[ip2*3+0], aXYZ[ip2*3+1], aXYZ[ip2*3+2]);
    const CVector3 p3(aXYZ[ip3*3+0], aXYZ[ip3*3+1], aXYZ[ip3*3+2]);
    double v0 = volume_Tet(p, p1,p2,p3);
    double v1 = volume_Tet(p0,p, p2,p3);
    double v2 = volume_Tet(p0,p1,p, p3);
    double v3 = volume_Tet(p0,p1,p2,p );
    double vt = (v0+v1+v2+v3);
    if (v0>-eps*vt && v1>-eps*vt && v2>-eps*vt && v3>-eps*vt ){
      double r0 = v0/(v0+v1+v2+v3);
      double r1 = v1/(v0+v1+v2+v3);
      double r2 = v2/(v0+v1+v2+v3);
      CPointElemSolid pt(itet1,r0,r1,r2);
      return pt;
    }
    if(      v0<v1 && v0<v2 && v0<v3 ){ itet1 = aTetSurRel[itet1*8+0*2+0]; }
    else if( v1<v0 && v1<v2 && v1<v3 ){ itet1 = aTetSurRel[itet1*8+1*2+0]; }
    else if( v2<v0 && v2<v1 && v2<v3 ){ itet1 = aTetSurRel[itet1*8+2*2+0]; }
    else{                               itet1 = aTetSurRel[itet1*8+3*2+0]; }
  }
  ////
  return CPointElemSolid();
}



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

CPointElemSurf nearest_Point_MeshTetFace3D
(const CVector3& p0,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTet,
 const std::vector<int>& aTetFaceSrf)
{
  const int noelTetFace[4][3] ={
    { 1, 2, 3 },
    { 0, 3, 2 },
    { 0, 1, 3 },
    { 0, 2, 1 } };
  ////
  double dist_min=-1.0;
  int itf_min = -1;
  CVector3 p_min;
  for(int itf=0;itf<aTetFaceSrf.size()/2;++itf){
    int itet = aTetFaceSrf[itf*2+0];
    int iface = aTetFaceSrf[itf*2+1];
    const int i0 = aTet[itet*4+noelTetFace[iface][0]];
    const int i1 = aTet[itet*4+noelTetFace[iface][1]];
    const int i2 = aTet[itet*4+noelTetFace[iface][2]];
    CVector3 q0 = CVector3(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2])-p0;
    CVector3 q1 = CVector3(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2])-p0;
    CVector3 q2 = CVector3(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2])-p0;
    double r0,r1;
    CVector3 p2 = nearest_Origin_Tri(r0,r1, q0,q1,q2);
    double dist = p2.Length();
    if( itf_min == -1 || dist < dist_min ){
      dist_min = dist;
      itf_min = itf;
      p_min = p2;
    }
  }
  assert( itf_min != -1 );
  {
    int itet = aTetFaceSrf[itf_min*2+0];
    int iface = aTetFaceSrf[itf_min*2+1];
    const int i0 = aTet[itet*4+noelTetFace[iface][0]];
    const int i1 = aTet[itet*4+noelTetFace[iface][1]];
    const int i2 = aTet[itet*4+noelTetFace[iface][2]];
    CVector3 q0(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
    CVector3 q1(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
    CVector3 q2(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
    double a0 = TriArea(p_min, q1,q2);
    double a1 = TriArea(p_min, q2,q0);
    double a2 = TriArea(p_min, q0,q1);
    int inva = 1.0/(a0+a1+a2);
    a0 *= inva;
    a1 *= inva;
    a2 *= inva;
    CPointElemSurf ptf;
    ptf.itri = itf_min;
    ptf.r0 = a0;
    ptf.r1 = a1;
    return ptf;
  }
}
