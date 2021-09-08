/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <map>

#include "delfem2/srchuni_v3.h"
#include "delfem2/geoproximity3_v3.h"
#include "delfem2/vec3.h"

namespace delfem2 {
namespace srchuni {

//! Volume of a tetrahedra
template<typename T>
T Volume_Tet(
    const CVec3<T> &v0,
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3) {
  const double v =
      (v1.p[0] - v0.p[0]) * ((v2.p[1] - v0.p[1]) * (v3.p[2] - v0.p[2]) - (v3.p[1] - v0.p[1]) * (v2.p[2] - v0.p[2]))
          + (v1.p[1] - v0.p[1])
              * ((v2.p[2] - v0.p[2]) * (v3.p[0] - v0.p[0]) - (v3.p[2] - v0.p[2]) * (v2.p[0] - v0.p[0]))
          + (v1.p[2] - v0.p[2])
              * ((v2.p[0] - v0.p[0]) * (v3.p[1] - v0.p[1]) - (v3.p[0] - v0.p[0]) * (v2.p[1] - v0.p[1]));
  return v * 0.16666666666666666666666666666667;
}

}
}

// ----------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::CPtElm3<T>::getPos_Tet
    (const std::vector<double> &aXYZ_,
     const std::vector<int> &aTet_) const {
  assert(ielem >= 0 && ielem < (int) aTet_.size() / 4);
  int ip0 = aTet_[ielem * 4 + 0];
  int ip1 = aTet_[ielem * 4 + 1];
  int ip2 = aTet_[ielem * 4 + 2];
  int ip3 = aTet_[ielem * 4 + 3];
  const CVec3<T> p0(aXYZ_[ip0 * 3 + 0], aXYZ_[ip0 * 3 + 1], aXYZ_[ip0 * 3 + 2]);
  const CVec3<T> p1(aXYZ_[ip1 * 3 + 0], aXYZ_[ip1 * 3 + 1], aXYZ_[ip1 * 3 + 2]);
  const CVec3<T> p2(aXYZ_[ip2 * 3 + 0], aXYZ_[ip2 * 3 + 1], aXYZ_[ip2 * 3 + 2]);
  const CVec3<T> p3(aXYZ_[ip3 * 3 + 0], aXYZ_[ip3 * 3 + 1], aXYZ_[ip3 * 3 + 2]);
  return r0 * p0 + r1 * p1 + r2 * p2 + (1.0 - r0 - r1 - r2) * p3;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::CPtElm3<double>::getPos_Tet
    (const std::vector<double> &aXYZ,
     const std::vector<int> &aTet) const;
#endif

// ----------------------------------------

template<typename T>
void delfem2::CPtElm3<T>::setPos_Tet(
    int it0,
    const CVec3<T> &q,
    const std::vector<double> &aXYZ_,
    const std::vector<int> &aTet_) {
  namespace lcl = delfem2::srchuni;
  assert(it0 >= 0 && it0 < (int) aTet_.size() / 4);
  int ip0 = aTet_[it0 * 4 + 0];
  int ip1 = aTet_[it0 * 4 + 1];
  int ip2 = aTet_[it0 * 4 + 2];
  int ip3 = aTet_[it0 * 4 + 3];
  const CVec3<T> p0(aXYZ_[ip0 * 3 + 0], aXYZ_[ip0 * 3 + 1], aXYZ_[ip0 * 3 + 2]);
  const CVec3<T> p1(aXYZ_[ip1 * 3 + 0], aXYZ_[ip1 * 3 + 1], aXYZ_[ip1 * 3 + 2]);
  const CVec3<T> p2(aXYZ_[ip2 * 3 + 0], aXYZ_[ip2 * 3 + 1], aXYZ_[ip2 * 3 + 2]);
  const CVec3<T> p3(aXYZ_[ip3 * 3 + 0], aXYZ_[ip3 * 3 + 1], aXYZ_[ip3 * 3 + 2]);
  double v0 = lcl::Volume_Tet(q, p1, p2, p3);
  double v1 = lcl::Volume_Tet(p0, q, p2, p3);
  double v2 = lcl::Volume_Tet(p0, p1, q, p3);
  //    double v3 = volume_Tet(p0,p1,p2, q);
  double vt = lcl::Volume_Tet(p0, p1, p2, p3);
  this->ielem = it0;
  this->r0 = v0 / vt;
  this->r1 = v1 / vt;
  this->r2 = v2 / vt;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CPtElm3<double>::setPos_Tet
    (int, const CVec3d &,
     const std::vector<double> &,
     const std::vector<int> &);
#endif

// --------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::CPtElm2<T>::Pos_Tri(
    const std::vector<double> &aXYZ_,
    const std::vector<unsigned int> &aTri_) const {
  assert(itri < aTri_.size() / 3);
  const int i0 = aTri_[itri * 3 + 0];
  const int i1 = aTri_[itri * 3 + 1];
  const int i2 = aTri_[itri * 3 + 2];
  const CVec3<T> p0(aXYZ_[i0 * 3 + 0], aXYZ_[i0 * 3 + 1], aXYZ_[i0 * 3 + 2]);
  const CVec3<T> p1(aXYZ_[i1 * 3 + 0], aXYZ_[i1 * 3 + 1], aXYZ_[i1 * 3 + 2]);
  const CVec3<T> p2(aXYZ_[i2 * 3 + 0], aXYZ_[i2 * 3 + 1], aXYZ_[i2 * 3 + 2]);
  return r0 * p0 + r1 * p1 + (1.0 - r0 - r1) * p2;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::CPtElm2<double>::Pos_Tri(
    const std::vector<double> &,
    const std::vector<unsigned int> &) const;
#endif

// --------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::CPtElm2<T>::Pos_Tri(
    const double *aXYZ_,
    [[maybe_unused]] unsigned int nXYZ,
    const unsigned int *aTri_,
    [[maybe_unused]] unsigned int nTri) const {
  assert(itri < nTri);
  const unsigned int i0 = aTri_[itri * 3 + 0];
  const unsigned int i1 = aTri_[itri * 3 + 1];
  const unsigned int i2 = aTri_[itri * 3 + 2];
  const CVec3<T> p0(aXYZ_[i0 * 3 + 0], aXYZ_[i0 * 3 + 1], aXYZ_[i0 * 3 + 2]);
  const CVec3<T> p1(aXYZ_[i1 * 3 + 0], aXYZ_[i1 * 3 + 1], aXYZ_[i1 * 3 + 2]);
  const CVec3<T> p2(aXYZ_[i2 * 3 + 0], aXYZ_[i2 * 3 + 1], aXYZ_[i2 * 3 + 2]);
  return r0 * p0 + r1 * p1 + (1.0 - r0 - r1) * p2;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::CPtElm2<double>::Pos_Tri(
    const double *,
    unsigned int,
    const unsigned int *,
    unsigned int) const;
#endif

// -----------------------------------------------


template<typename T>
delfem2::CVec3<T> delfem2::CPtElm2<T>::UNorm_Tri(
	[[maybe_unused]] const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri0,
    const std::vector<double> &aNorm0) const {
  assert(itri < aTri0.size() / 3);
  const unsigned int i0 = aTri0[itri * 3 + 0];
  const unsigned int i1 = aTri0[itri * 3 + 1];
  const unsigned int i2 = aTri0[itri * 3 + 2];
  const CVec3<T> n0(aNorm0[i0 * 3 + 0], aNorm0[i0 * 3 + 1], aNorm0[i0 * 3 + 2]);
  const CVec3<T> n1(aNorm0[i1 * 3 + 0], aNorm0[i1 * 3 + 1], aNorm0[i1 * 3 + 2]);
  const CVec3<T> n2(aNorm0[i2 * 3 + 0], aNorm0[i2 * 3 + 1], aNorm0[i2 * 3 + 2]);
  return (r0 * n0 + r1 * n1 + (1.0 - r0 - r1) * n2).normalized();
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::CPtElm2<double>::UNorm_Tri(
    const std::vector<double> &,
    const std::vector<unsigned int> &,
    const std::vector<double> &) const;
#endif

// ------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::CPtElm2<T>::UNorm_Tri(
    [[maybe_unused]] const double *aXYZ,
    [[maybe_unused]] unsigned int nXYZ,
    const unsigned int *aTri,
    [[maybe_unused]] unsigned int nTri,
    const double *aNorm) const {
  assert(itri < nTri);
  const unsigned int i0 = aTri[itri * 3 + 0];
  const unsigned int i1 = aTri[itri * 3 + 1];
  const unsigned int i2 = aTri[itri * 3 + 2];
  const CVec3<T> n0(aNorm[i0 * 3 + 0], aNorm[i0 * 3 + 1], aNorm[i0 * 3 + 2]);
  const CVec3<T> n1(aNorm[i1 * 3 + 0], aNorm[i1 * 3 + 1], aNorm[i1 * 3 + 2]);
  const CVec3<T> n2(aNorm[i2 * 3 + 0], aNorm[i2 * 3 + 1], aNorm[i2 * 3 + 2]);
  return (r0 * n0 + r1 * n1 + (1.0 - r0 - r1) * n2).normalized();
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::CPtElm2<double>::UNorm_Tri(
    const double *aXYZ,
    unsigned int nXYZ,
    const unsigned int *aTri,
    unsigned int nTri,
    const double *aNorm) const;
#endif

// ------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::CPtElm2<T>::Pos_TetFace(
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTet,
    const std::vector<int> &aTetFace) const {
  const int noelTetFace[4][3] = {
      {1, 2, 3},
      {0, 3, 2},
      {0, 1, 3},
      {0, 2, 1}};
  int itet = aTetFace[itri * 2 + 0];
  int iface = aTetFace[itri * 2 + 1];
  double r2 = (1 - r0 - r1);
  int ielno0 = noelTetFace[iface][0];
  int ielno1 = noelTetFace[iface][1];
  int ielno2 = noelTetFace[iface][2];
  int iq0 = aTet[itet * 4 + ielno0];
  int iq1 = aTet[itet * 4 + ielno1];
  int iq2 = aTet[itet * 4 + ielno2];
  CVec3<T> p;
  p.p[0] = r0 * aXYZ[iq0 * 3 + 0] + r1 * aXYZ[iq1 * 3 + 0] + r2 * aXYZ[iq2 * 3 + 0];
  p.p[1] = r0 * aXYZ[iq0 * 3 + 1] + r1 * aXYZ[iq1 * 3 + 1] + r2 * aXYZ[iq2 * 3 + 1];
  p.p[2] = r0 * aXYZ[iq0 * 3 + 2] + r1 * aXYZ[iq1 * 3 + 2] + r2 * aXYZ[iq2 * 3 + 2];
  return p;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::CPtElm2<double>::Pos_TetFace(
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTet,
    const std::vector<int> &aTetFace) const;
#endif

// ----------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::CPtElm2<T>::Pos_Grid(
    unsigned int nx,
    unsigned int ny,
    double el,
    std::vector<float> &aH) const {
  int iey = (this->itri / 2) / nx;
  int iex = (this->itri / 2) - nx * iey;
  CVec3<T> p00((iex + 0) * el, (iey + 0) * el, aH[(iey + 0) * nx + (iex + 0)]);
  CVec3<T> p10((iex + 1) * el, (iey + 0) * el, aH[(iey + 0) * nx + (iex + 1)]);
  CVec3<T> p01((iex + 0) * el, (iey + 1) * el, aH[(iey + 1) * nx + (iex + 0)]);
  CVec3<T> p11((iex + 1) * el, (iey + 1) * el, aH[(iey + 1) * nx + (iex + 1)]);
  if (this->itri % 2 == 0) { return p00 * r0 + p10 * r1 + p11 * (1 - r0 - r1); }
  return p00 * r0 + p11 * r1 + p01 * (1 - r0 - r1);
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::CPtElm2<double>::Pos_Grid(
    unsigned int nx, unsigned int ny,
    double el, std::vector<float> &aH) const;
#endif

// ----------------------------------------

template<typename T>
bool delfem2::CPtElm2<T>::Check(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    double eps) const {
  if (itri >= aTri.size() / 3) { return false; }
  if (r0 < -eps || r0 > 1 + eps) { return false; }
  if (r1 < -eps || r1 > 1 + eps) { return false; }
  double r2 = 1 - r0 - r1;
  if (r2 < -eps || r2 > 1 + eps) { return false; }
  return true;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::CPtElm2<double>::Check(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    double eps) const;
#endif

// ------------------------------------------



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

// ---------------------------------------------------------------------------------

template<typename T>
std::vector<delfem2::CPtElm2<T>>
delfem2::IntersectionLine_MeshTri3(
    const CVec3<T> &org,
    const CVec3<T> &dir,
    const std::vector<unsigned int> &aTri_,
    const std::vector<T> &aXYZ_,
    T eps) {
  std::vector<CPtElm2<T>> aPES;
  for (unsigned int itri = 0; itri < aTri_.size() / 3; ++itri) {
    const unsigned int ip0 = aTri_[itri * 3 + 0];
    const unsigned int ip1 = aTri_[itri * 3 + 1];
    const unsigned int ip2 = aTri_[itri * 3 + 2];
    assert(ip0 < aXYZ_.size() / 3);
    assert(ip1 < aXYZ_.size() / 3);
    assert(ip2 < aXYZ_.size() / 3);
    const CVec3<T> p0(aXYZ_[ip0 * 3 + 0], aXYZ_[ip0 * 3 + 1], aXYZ_[ip0 * 3 + 2]);
    const CVec3<T> p1(aXYZ_[ip1 * 3 + 0], aXYZ_[ip1 * 3 + 1], aXYZ_[ip1 * 3 + 2]);
    const CVec3<T> p2(aXYZ_[ip2 * 3 + 0], aXYZ_[ip2 * 3 + 1], aXYZ_[ip2 * 3 + 2]);
    double r0, r1;
    bool res = delfem2::IntersectRay_Tri3(r0, r1,
                                          org, dir, p0, p1, p2,
                                          eps);
    if (!res) { continue; }
    aPES.emplace_back(itri, r0, r1);
  }
  return aPES;
}
#ifdef DFM2_STATIC_LIBRARY
template
std::vector<delfem2::CPtElm2<double>> delfem2::IntersectionLine_MeshTri3(
    const CVec3<double> &org,
    const CVec3<double> &dir,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ,
    double eps);
#endif

// -------------------------------------

template<typename T>
void delfem2::IntersectionRay_MeshTri3(
    std::map<T, CPtElm2<T>> &mapDepthPES,
    const CVec3<T> &org,
    const CVec3<T> &dir,
    const std::vector<unsigned int> &aTri_,
    const std::vector<T> &aXYZ_,
    T eps) {
  const std::vector<CPtElm2<T>> aPES = IntersectionLine_MeshTri3(
      org, dir,
      aTri_, aXYZ_,
      eps);
  mapDepthPES.clear();
  for (auto pes : aPES) {
    CVec3<T> p0 = pes.Pos_Tri(aXYZ_, aTri_);
    double depth = (p0 - org).dot(dir);
    if (depth < 0) continue;
    mapDepthPES.insert(std::make_pair(depth, pes));
  }
}
#ifdef DFM2_STATIC_LIBRARY
template
void delfem2::IntersectionRay_MeshTri3(
    std::map<double, CPtElm2<double>> &mapDepthPES,
    const CVec3<double> &org, const CVec3<double> &dir,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ,
    double eps);
#endif

// ----------------------

template<typename T>
void delfem2::IntersectionRay_MeshTri3DPart(
    std::map<T, CPtElm2<T>> &mapDepthPES,
    const CVec3<T> &org,
    const CVec3<T> &dir,
    const std::vector<unsigned int> &aTri_,
    const std::vector<T> &aXYZ_,
    const std::vector<unsigned int> &aIndTri,
    T eps) {
  mapDepthPES.clear();
  for (unsigned int itri : aIndTri) {
    const unsigned int ip0 = aTri_[itri * 3 + 0];
    assert(ip0 < aXYZ_.size() / 3);
    const unsigned int ip1 = aTri_[itri * 3 + 1];
    assert(ip1 < aXYZ_.size() / 3);
    const unsigned int ip2 = aTri_[itri * 3 + 2];
    assert(ip2 < aXYZ_.size() / 3);
    const CVec3<T> p0(aXYZ_[ip0 * 3 + 0], aXYZ_[ip0 * 3 + 1], aXYZ_[ip0 * 3 + 2]);
    const CVec3<T> p1(aXYZ_[ip1 * 3 + 0], aXYZ_[ip1 * 3 + 1], aXYZ_[ip1 * 3 + 2]);
    const CVec3<T> p2(aXYZ_[ip2 * 3 + 0], aXYZ_[ip2 * 3 + 1], aXYZ_[ip2 * 3 + 2]);
    double r0, r1;
    bool res = IntersectRay_Tri3(r0, r1,
                                 org, dir, p0, p1, p2,
                                 eps);
    if (!res) { continue; }
    double r2 = 1 - r0 - r1;
    CVec3<T> q0 = p0 * r0 + p1 * r1 + p2 * r2;
    double depth = (q0 - org).dot(dir) / dir.squaredNorm();
    if (depth < 0) continue;
    mapDepthPES.insert(std::make_pair(depth, CPtElm2<T>(itri, r0, r1)));
  }
}
#ifdef DFM2_STATIC_LIBRARY
template
void delfem2::IntersectionRay_MeshTri3DPart(
    std::map<double, CPtElm2<double>> &,
    const CVec3<double> &,
    const CVec3<double> &,
    const std::vector<unsigned int> &,
    const std::vector<double> &,
    const std::vector<unsigned int> &,
    double);
#endif

// -----------------------------------------------

/*
CPtElm2 intersect_Ray_MeshTriFlag3D
(const CVector3& org, const CVector3& dir,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ,
 int iflag,
 const std::vector<int>& aFlag)
{
  assert( aTri.size()/3 == aFlag.size() );
  std::map<double, CPtElm2> pickMap;
  for (int itri = 0; itri<(int)aTri.size()/3; itri++){
    if( aFlag[itri] != iflag ) continue;
    double depth;
    CPtElm2 res = intersect_Ray_Tri3D(depth, org,dir, itri, aTri,aXYZ);
    if( res.itri<0 ){ continue; }
    pickMap.insert(std::make_pair(depth, res));
  }
  if (pickMap.empty()) return CPtElm2();
  return pickMap.begin()->second;
}

CPtElm2 intersect_Ray_MeshTri3D
(const CVector3& org, const CVector3& dir,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ)
{
  std::map<double, CPtElm2> pickMap;
  for (int itri = 0; itri<(int)aTri.size()/3; itri++){
    double depth;
    CPtElm2 res = intersect_Ray_Tri3D(depth, org,dir, itri, aTri,aXYZ);
    if( res.itri<0 ){ continue; }
    pickMap.insert(std::make_pair(depth, res));
  }
  if (pickMap.empty()) return CPtElm2(-1,0,0);
  return pickMap.begin()->second;
}

CPtElm2 intersect_Ray_MeshTri3D
(const CVector3& org, const CVector3& dir,
 int itri_start, // starting triangle
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTriSurRel)
{
  int itri1 = itri_start;
  if (itri1<0||itri1>=(int)aTri.size()/3){ return CPtElm2(); }
  for (int itr = 0; itr<50; itr++){
    if (itri1==-1) return CPtElm2();
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
      return CPtElm2(itri1,r0,r1);
    }
    if (v0<v1 && v0<v2){      itri1 = aTriSurRel[itri1*6+0*2+0]; }
    else if (v1<v0 && v1<v2){ itri1 = aTriSurRel[itri1*6+1*2+0]; }
    else{                     itri1 = aTriSurRel[itri1*6+2*2+0]; }
  }
  ////
  return CPtElm2();
}

CPtElm2 intersect_Ray_MeshTri3D
(const CVector3& org, const CVector3& dir,
 const CPtElm2& ptri0,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTriSurRel)
{
  CPtElm2 ptri;
  if( ptri0.itri != -1 ){
    ptri = intersect_Ray_MeshTri3D(org,dir, ptri0.itri,aTri,aXYZ,aTriSurRel);
  }
  if( ptri.itri == -1 ){
    ptri = intersect_Ray_MeshTri3D(org,dir,aTri,aXYZ);
  }
  if( ptri.itri == -1 ) return ptri;
  return ptri;
}
 */


DFM2_INLINE void delfem2::IntersectionLine_Hightfield(
    std::vector<CPtElm2<double>> &aPes,
    const double src[3],
    const double dir[3],
    unsigned int nx,
    unsigned int ny,
    const std::vector<float> &aH) {
  assert(aH.size() == nx * ny);
  for (unsigned int iey = 0; iey < ny - 1; ++iey) {
    for (unsigned int iex = 0; iex < nx - 1; ++iex) {
      const double h00 = aH[(iey + 0) * nx + (iex + 0)];
      const double h10 = aH[(iey + 0) * nx + (iex + 1)];
      const double h01 = aH[(iey + 1) * nx + (iex + 0)];
      const double h11 = aH[(iey + 1) * nx + (iex + 1)];
      const double p00[3] = {double(iex + 0), double(iey + 0), h00};
      const double p10[3] = {double(iex + 1), double(iey + 0), h10};
      const double p01[3] = {double(iex + 0), double(iey + 1), h01};
      const double p11[3] = {double(iex + 1), double(iey + 1), h11};
      double r0 = 1.0, r1 = 1.0;
      if (IntersectRay_Tri3(r0, r1,
                            CVec3d(src), CVec3d(dir),
                            CVec3d(p00), CVec3d(p10), CVec3d(p11), 1.0e-3)) {
        aPes.emplace_back((iey * nx + iex) * 2 + 0, r0, r1);
      }
      // ---------------------
      if (IntersectRay_Tri3(r0, r1,
                            CVec3d(src), CVec3d(dir),
                            CVec3d(p00), CVec3d(p11), CVec3d(p01), 1.0e-3)) {
        aPes.emplace_back((iey * nx + iex) * 2 + 1, r0, r1);
      }
    } // iex
  } // iey
}

// ----------------------------------------------------------

template<typename T>
delfem2::CPtElm2<T> delfem2::Nearest_Point_MeshTri3D(
    const CVec3<T> &q,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri) {
  CPtElm2<T> pes;
  double min_dist = -1;
  const size_t nTri = aTri.size() / 3;
  for (unsigned int it = 0; it < nTri; ++it) {
    const unsigned int i0 = aTri[it * 3 + 0];
    const unsigned int i1 = aTri[it * 3 + 1];
    const unsigned int i2 = aTri[it * 3 + 2];
    const CVec3<T> p0(aXYZ[i0 * 3 + 0] - q.x, aXYZ[i0 * 3 + 1] - q.y, aXYZ[i0 * 3 + 2] - q.z);
    const CVec3<T> p1(aXYZ[i1 * 3 + 0] - q.x, aXYZ[i1 * 3 + 1] - q.y, aXYZ[i1 * 3 + 2] - q.z);
    const CVec3<T> p2(aXYZ[i2 * 3 + 0] - q.x, aXYZ[i2 * 3 + 1] - q.y, aXYZ[i2 * 3 + 2] - q.z);
    double r0, r1;
    CVec3<T> p_min = Nearest_Origin_Tri(r0, r1, p0, p1, p2);
    double dist = p_min.squaredNorm();
    if (min_dist < 0 || dist < min_dist) {
      min_dist = dist;
      pes = CPtElm2<T>((int) it, r0, r1);
    }
  }
  assert(pes.itri != UINT_MAX);
  return pes;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CPtElm2<double> delfem2::Nearest_Point_MeshTri3D(
    const CVec3d &,
    const std::vector<double> &,
    const std::vector<unsigned int> &);
#endif

// ----------------------------------------------


template<typename T>
delfem2::CPtElm2<T> delfem2::Nearest_Point_MeshTri3DPart(
    const CVec3<T> &q,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<int> &aIndTri_Cand) {
  double min_dist = -1;
  CPtElm2<T> pes;
  for (int itri0 : aIndTri_Cand) {
    const unsigned int i0 = aTri[itri0 * 3 + 0];
    const unsigned int i1 = aTri[itri0 * 3 + 1];
    const unsigned int i2 = aTri[itri0 * 3 + 2];
    const CVec3<T> p0(aXYZ[i0 * 3 + 0] - q.x, aXYZ[i0 * 3 + 1] - q.y, aXYZ[i0 * 3 + 2] - q.z);
    const CVec3<T> p1(aXYZ[i1 * 3 + 0] - q.x, aXYZ[i1 * 3 + 1] - q.y, aXYZ[i1 * 3 + 2] - q.z);
    const CVec3<T> p2(aXYZ[i2 * 3 + 0] - q.x, aXYZ[i2 * 3 + 1] - q.y, aXYZ[i2 * 3 + 2] - q.z);
    double r0, r1;
    CVec3<T> p_min = Nearest_Origin_Tri(r0, r1, p0, p1, p2);
    assert(r0 > -1.0e-10 && r1 > -1.0e-10 && (1 - r0 - r1) > -1.0e-10);
    double dist = p_min.squaredNorm();
    if (min_dist < 0 || dist < min_dist) {
      min_dist = dist;
      pes = CPtElm2<T>(itri0, r0, r1);
    }
  }
  return pes;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CPtElm2<double> delfem2::Nearest_Point_MeshTri3DPart(
    const CVec3d &,
    const std::vector<double> &,
    const std::vector<unsigned int> &,
    const std::vector<int> &);
#endif

// ----------------------------------------------------------------------------

template<typename T>
delfem2::CPtElm3<T> delfem2::Nearest_Point_MeshTet3D
    (const CVec3<T> &q,
     const std::vector<double> &aXYZ,
     const std::vector<int> &aTet) {
  const double eps = 1.0e-4;
  const unsigned int ntet = aTet.size() / 4;
  for (unsigned int itet = 0; itet < ntet; ++itet) {
    CPtElm3<T> pt;
    pt.setPos_Tet(itet, q, aXYZ, aTet);
    if (pt.isInside(-eps)) { return pt; }
  }
  return CPtElm3<T>();
}

// ------------------------------

template<typename T>
delfem2::CPtElm3<T> delfem2::Nearest_Point_MeshTet3D
    (const CVec3<T> &p,
     int itet_start, // starting triangle
     const std::vector<double> &aXYZ,
     const std::vector<int> &aTet,
     const std::vector<int> &aTetSurRel) {
  const double eps = 1.0e-4;
  int itet1 = itet_start;
  if (itet1 < 0 || itet1 >= (int) aTet.size() / 4) { return CPtElm3<T>(); }
  for (int itr = 0; itr < 50; itr++) {
    if (itet1 == -1) return CPtElm3<T>();
    int ip0 = aTet[itet1 * 4 + 0];
    int ip1 = aTet[itet1 * 4 + 1];
    int ip2 = aTet[itet1 * 4 + 2];
    int ip3 = aTet[itet1 * 4 + 3];
    const CVec3<T> p0(aXYZ[ip0 * 3 + 0], aXYZ[ip0 * 3 + 1], aXYZ[ip0 * 3 + 2]);
    const CVec3<T> p1(aXYZ[ip1 * 3 + 0], aXYZ[ip1 * 3 + 1], aXYZ[ip1 * 3 + 2]);
    const CVec3<T> p2(aXYZ[ip2 * 3 + 0], aXYZ[ip2 * 3 + 1], aXYZ[ip2 * 3 + 2]);
    const CVec3<T> p3(aXYZ[ip3 * 3 + 0], aXYZ[ip3 * 3 + 1], aXYZ[ip3 * 3 + 2]);
    double v0 = Volume_Tet(p, p1, p2, p3);
    double v1 = Volume_Tet(p0, p, p2, p3);
    double v2 = Volume_Tet(p0, p1, p, p3);
    double v3 = Volume_Tet(p0, p1, p2, p);
    double vt = (v0 + v1 + v2 + v3);
    if (v0 > -eps * vt && v1 > -eps * vt && v2 > -eps * vt && v3 > -eps * vt) {
      double r0 = v0 / (v0 + v1 + v2 + v3);
      double r1 = v1 / (v0 + v1 + v2 + v3);
      double r2 = v2 / (v0 + v1 + v2 + v3);
      CPtElm3<T> pt(itet1, r0, r1, r2);
      return pt;
    }
    if (v0 < v1 && v0 < v2 && v0 < v3) { itet1 = aTetSurRel[itet1 * 8 + 0 * 2 + 0]; }
    else if (v1 < v0 && v1 < v2 && v1 < v3) { itet1 = aTetSurRel[itet1 * 8 + 1 * 2 + 0]; }
    else if (v2 < v0 && v2 < v1 && v2 < v3) { itet1 = aTetSurRel[itet1 * 8 + 2 * 2 + 0]; }
    else { itet1 = aTetSurRel[itet1 * 8 + 3 * 2 + 0]; }
  }
  return CPtElm3<T>();
}

// ---------------------------------------------------------------

template<typename T>
delfem2::CPtElm2<T> delfem2::Nearest_Point_MeshTetFace3D
    (const CVec3<T> &p0,
     const std::vector<double> &aXYZ,
     const std::vector<int> &aTet,
     const std::vector<int> &aTetFaceSrf) {
  const int noelTetFace[4][3] = {
      {1, 2, 3},
      {0, 3, 2},
      {0, 1, 3},
      {0, 2, 1}};
  ////
  double dist_min = -1.0;
  int itf_min = -1;
  CVec3<T> p_min;
  for (size_t itf = 0; itf < aTetFaceSrf.size() / 2; ++itf) {
    int itet = aTetFaceSrf[itf * 2 + 0];
    int iface = aTetFaceSrf[itf * 2 + 1];
    const int i0 = aTet[itet * 4 + noelTetFace[iface][0]];
    const int i1 = aTet[itet * 4 + noelTetFace[iface][1]];
    const int i2 = aTet[itet * 4 + noelTetFace[iface][2]];
    CVec3<T> q0 = CVec3<T>(aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]) - p0;
    CVec3<T> q1 = CVec3<T>(aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]) - p0;
    CVec3<T> q2 = CVec3<T>(aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]) - p0;
    double r0, r1;
    CVec3<T> p2 = Nearest_Origin_Tri(r0, r1, q0, q1, q2);
    double dist = p2.norm();
    if (itf_min == -1 || dist < dist_min) {
      dist_min = dist;
      itf_min = itf;
      p_min = p2;
    }
  }
  assert(itf_min != -1);
  {
    int itet = aTetFaceSrf[itf_min * 2 + 0];
    int iface = aTetFaceSrf[itf_min * 2 + 1];
    const int i0 = aTet[itet * 4 + noelTetFace[iface][0]];
    const int i1 = aTet[itet * 4 + noelTetFace[iface][1]];
    const int i2 = aTet[itet * 4 + noelTetFace[iface][2]];
    CVec3<T> q0(aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]);
    CVec3<T> q1(aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]);
    CVec3<T> q2(aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]);
    double a0 = Area_Tri(p_min, q1, q2);
    double a1 = Area_Tri(p_min, q2, q0);
    double a2 = Area_Tri(p_min, q0, q1);
    double inva = 1.0 / (a0 + a1 + a2);
    a0 *= inva;
    a1 *= inva;
    a2 *= inva;
    CPtElm2<T> ptf;
    ptf.itri = itf_min;
    ptf.r0 = a0;
    ptf.r1 = a1;
    return ptf;
  }
}

// ----------------------------------------

template<typename T>
double delfem2::SDFNormal_NearestPoint
    (CVec3<T> &n0,
     const CVec3<T> &p0,
     const CPtElm2<T> &pes,
     const double *aXYZ, unsigned int nXYZ,
     const unsigned int *aTri, unsigned int nTri,
     const double *aNorm) {
  CVec3<T> q1 = pes.Pos_Tri(aXYZ, nXYZ, aTri, nTri);
  double dist = (q1 - p0).norm();
  CVec3<T> n1 = pes.UNorm_Tri(aXYZ, nXYZ, aTri, nTri, aNorm);
  if ((q1 - p0).dot(n1) > 0) {  //inside
    if (dist < 1.0e-6) { n0 = n1; }
    else { n0 = (q1 - p0).normalized(); }
    return dist;
  } else { // outside
    if (dist < 1.0e-6) { n0 = n1; }
    else { n0 = (p0 - q1).normalized(); }
    return -dist;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::SDFNormal_NearestPoint
    (CVec3d &n0,
     const CVec3d &p0,
     const CPtElm2<double> &pes,
     const double *aXYZ, unsigned int nXYZ,
     const unsigned int *aTri, unsigned int nTri,
     const double *aNorm);
#endif

// ------------------------------------------------------

template<typename T>
double delfem2::SDFNormal_NearestPoint
    (CVec3<T> &n0,
     const CVec3<T> &p0,
     const CPtElm2<T> &pes,
     const std::vector<double> &aXYZ,
     const std::vector<unsigned int> &aTri,
     const std::vector<double> &aNorm) {
  CVec3<T> q1 = pes.Pos_Tri(aXYZ, aTri);
  double dist = (q1 - p0).norm();
  CVec3<T> n1 = pes.UNorm_Tri(aXYZ, aTri, aNorm);
  if ((q1 - p0).dot(n1) > 0) {  //inside
    if (dist < 1.0e-6) { n0 = n1; }
    else { n0 = (q1 - p0).normalized(); }
    return dist;
  } else { // outside
    if (dist < 1.0e-6) { n0 = n1; }
    else { n0 = (p0 - q1).normalized(); }
    return -dist;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::SDFNormal_NearestPoint
    (CVec3d &n0,
     const CVec3d &p0,
     const CPtElm2<double> &pes,
     const std::vector<double> &aXYZ,
     const std::vector<unsigned int> &aTri,
     const std::vector<double> &aNorm);
#endif

// ---------------------------------------------

template<typename T>
double delfem2::DistanceToTri
    (CPtElm2<T> &pes,
     const CVec3<T> &p,
     unsigned int itri0,
     const std::vector<double> &aXYZ_,
     const std::vector<unsigned int> &aTri_) {
  const unsigned int i0 = aTri_[itri0 * 3 + 0];
  const unsigned int i1 = aTri_[itri0 * 3 + 1];
  const unsigned int i2 = aTri_[itri0 * 3 + 2];
  const CVec3<T> p0(aXYZ_[i0 * 3 + 0] - p.p[0], aXYZ_[i0 * 3 + 1] - p.p[1], aXYZ_[i0 * 3 + 2] - p.p[2]);
  const CVec3<T> p1(aXYZ_[i1 * 3 + 0] - p.p[0], aXYZ_[i1 * 3 + 1] - p.p[1], aXYZ_[i1 * 3 + 2] - p.p[2]);
  const CVec3<T> p2(aXYZ_[i2 * 3 + 0] - p.p[0], aXYZ_[i2 * 3 + 1] - p.p[1], aXYZ_[i2 * 3 + 2] - p.p[2]);
  double r0, r1;
  CVec3<T> p_min = Nearest_Origin_Tri(r0, r1, p0, p1, p2);
  assert(r0 > -1.0e-10 && r1 > -1.0e-10 && (1 - r0 - r1) > -1.0e-10);
  pes.itri = itri0;
  pes.r0 = r0;
  pes.r1 = r1;
  return p_min.norm();
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::DistanceToTri
    (CPtElm2<double> &pes,
     const CVec3<double> &p,
     unsigned int itri0,
     const std::vector<double> &aXYZ,
     const std::vector<unsigned int> &aTri);
#endif

// -------------------------------------------

template<typename T>
double delfem2::DistanceToTri(
    CPtElm2<T> &pes,
    const CVec3<T> &p,
    unsigned int itri0,
    const double *aXYZ_,
    size_t nXYZ,
    const unsigned int *aTri_,
    size_t nTri) {
  const unsigned int i0 = aTri_[itri0 * 3 + 0];
  const unsigned int i1 = aTri_[itri0 * 3 + 1];
  const unsigned int i2 = aTri_[itri0 * 3 + 2];
  const CVec3<T> p0(aXYZ_[i0 * 3 + 0] - p.p[0], aXYZ_[i0 * 3 + 1] - p.p[1], aXYZ_[i0 * 3 + 2] - p.p[2]);
  const CVec3<T> p1(aXYZ_[i1 * 3 + 0] - p.p[0], aXYZ_[i1 * 3 + 1] - p.p[1], aXYZ_[i1 * 3 + 2] - p.p[2]);
  const CVec3<T> p2(aXYZ_[i2 * 3 + 0] - p.p[0], aXYZ_[i2 * 3 + 1] - p.p[1], aXYZ_[i2 * 3 + 2] - p.p[2]);
  double r0, r1;
  CVec3<T> p_min = Nearest_Origin_Tri(r0, r1, p0, p1, p2);
  assert(r0 > -1.0e-10 && r1 > -1.0e-10 && (1 - r0 - r1) > -1.0e-10);
  pes.itri = itri0;
  pes.r0 = r0;
  pes.r1 = r1;
  return p_min.norm();
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::DistanceToTri(
    CPtElm2<double> &pes,
    const CVec3<double> &p,
    unsigned int itri0,
    const double *aXYZ,
    size_t nXYZ,
    const unsigned int *aTri,
    size_t nTri);
#endif

