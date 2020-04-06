/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_SRCHUNI_V3_H
#define DFM2_SRCHUNI_V3_H

#include "delfem2/dfm2_inline.h"
#include <stdio.h>
#include <map>
#include "delfem2/vec3.h"

// -----------------------------------

namespace delfem2 {
  
template <typename T>
class CPointElemSolid{
public:
  CPointElemSolid(): ielem(-1), r0(0), r1(0), r2(0) {}
  CPointElemSolid(int itet, double r0, double r1, double r2):ielem(itet),r0(r0),r1(r1),r2(r2) {}
  bool isInside(double eps) const {
    double r3 = (1-r0-r1-r2);
    if( r0 > eps && r1 > eps && r2 > eps &&  r3 > eps ){ return true; }
    return false;
  }
  int indFace(double eps) const {
    double r3 = (1-r0-r1-r2);
    if( fabs(r0)<eps && fabs(r1)>eps && fabs(r2)>eps && fabs(r3)>eps ){ return 0; }
    if( fabs(r0)>eps && fabs(r1)<eps && fabs(r2)>eps && fabs(r3)>eps ){ return 1; }
    if( fabs(r0)>eps && fabs(r1)>eps && fabs(r2)<eps && fabs(r3)>eps ){ return 2; }
    if( fabs(r0)>eps && fabs(r1)>eps && fabs(r2)>eps && fabs(r3)<eps ){ return 3; }
    return -1;
  }
  delfem2::CVec3<T> getPos_Tet(const std::vector<double>& aXYZ, const std::vector<int>& aTet) const;
  void setPos_Tet(int it0, const delfem2::CVec3<T>& q, const std::vector<double>& aXYZ, const std::vector<int>& aTet);
public:
  int ielem;
  double r0, r1, r2;
};

// --------------------------------------------------

// class for point in a triangle
template <typename T>
class CPointElemSurf{
public:
  CPointElemSurf() : itri(-1), r0(0), r1(0) {}
  CPointElemSurf(int itri, double r0, double r1):itri(itri), r0(r0),r1(r1) {}
  delfem2::CVec3<T> Pos_Tri(const std::vector<double>& aXYZ,
                   const std::vector<unsigned int>& aTri) const;
  delfem2::CVec3<T> Pos_Tri(const double* aXYZ, unsigned int nXYZ,
                   const unsigned int* aTri, unsigned int nTri) const;
  delfem2::CVec3<T> Pos_TetFace(const std::vector<double>& aXYZ,
                       const std::vector<int>& aTet,
                       const std::vector<int>& aTetFace) const;
  delfem2::CVec3<T> Pos_Grid(
      unsigned int nx, unsigned int ny,
      double el,
      std::vector<float>& aH) const;
  delfem2::CVec3<T> UNorm_Tri(const std::vector<double>& aXYZ,
                     const std::vector<unsigned int>& aTri,
                     const std::vector<double>& aNorm) const;
  delfem2::CVec3<T> UNorm_Tri(const double* aXYZ, unsigned int nXYZ,
                     const unsigned int* aTri, unsigned int nTri,
                     const double* aNorm) const;
  bool Check(const std::vector<double>& aXYZ,
             const std::vector<unsigned int>& aTri,
             double eps) const;
public:
  int itri; // can be -1
  double r0, r1;
};
  
template <typename T>
std::ostream &operator << (std::ostream &output,
                           const CPointElemSurf<T>& v)
{
  output.setf(std::ios::scientific);
  output << v.itri << " " << v.r0 << " " << v.r1;
  return output;
}

template <typename T>
std::istream &operator >> (std::istream &input,
                           CPointElemSurf<T>& v){
  input>>v.itri>>v.r0>>v.r1;
  return input;
}

// ----------------------------------------------------------

template <typename REAL>
std::vector<CPointElemSurf<REAL>> IntersectionLine_MeshTri3(
    const delfem2::CVec3<REAL>& org, const delfem2::CVec3<REAL>& dir,
    const std::vector<unsigned int>& aTri,
    const std::vector<REAL>& aXYZ,
    REAL eps);


template <typename REAL>
void IntersectionRay_MeshTri3 (
    std::map<REAL,CPointElemSurf<REAL>>& mapDepthPES,
    const delfem2::CVec3<REAL>& org,
    const delfem2::CVec3<REAL>& dir,
    const std::vector<unsigned int>& aTri,
    const std::vector<REAL>& aXYZ,
    REAL eps );
  
template <typename REAL>
void IntersectionRay_MeshTri3DPart(
    std::map<REAL,CPointElemSurf<REAL>>& mapDepthPES,
    const delfem2::CVec3<REAL>& org, const delfem2::CVec3<REAL>& dir,
    const std::vector<unsigned int>& aTri,
    const std::vector<REAL>& aXYZ,
    const std::vector<int>& aIndTri,
    REAL eps);

/*
 CPointElemSurf intersect_Ray_Tri3D(double& depth,
 const delfem2::CVector3& org, const delfem2::CVector3& dir,
 int itri,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ);
 */
/*
CPointElemSurf intersect_Ray_MeshTri3D(const delfem2::CVector3& org, const delfem2::CVector3& dir,
                                       const std::vector<unsigned int>& aTri,
                                       const std::vector<double>& aXYZ);
CPointElemSurf intersect_Ray_MeshTri3D(const delfem2::CVector3& org, const delfem2::CVector3& dir,
                                       int itri_start, // starting triangle
                                       const std::vector<unsigned int>& aTri,
                                       const std::vector<double>& aXYZ,
                                       const std::vector<int>& aTriSurRel);
CPointElemSurf intersect_Ray_MeshTri3D(const delfem2::CVector3& org, const delfem2::CVector3& dir,
                                       const CPointElemSurf& ptri0,
                                       const std::vector<unsigned int>& aTri,
                                       const std::vector<double>& aXYZ,
                                       const std::vector<int>& aTriSurRel);
CPointElemSurf intersect_Ray_MeshTriFlag3D(const delfem2::CVector3& org, const delfem2::CVector3& dir,
                                           const std::vector<unsigned int>& aTri,
                                           const std::vector<double>& aXYZ,
                                           int iflag,
                                           const std::vector<int>& aFlag);
 */

DFM2_INLINE void IntersectionLine_Hightfield(
    std::vector<CPointElemSurf<double>>& aPos,
    //
    double hmin, double hmax,
    const double src[3],
    const double dir[3],
    double nx, double ny, double elen,
    const std::vector<float>& aH);

// above functions for ray interesection
// -----------------------------------------------------------
// below functions for nearest


template <typename T>
CPointElemSurf<T> Nearest_Point_MeshTri3D(const delfem2::CVec3<T>& q,
                                       const std::vector<double>& aXYZ,
                                       const std::vector<unsigned int>& aTri);

template <typename T>
CPointElemSurf<T> Nearest_Point_MeshTri3DPart(const delfem2::CVec3<T>& q,
                                           const std::vector<double>& aXYZ,
                                           const std::vector<unsigned int>& aTri,
                                           const std::vector<int>& aIndTri_Cand);

template <typename T>
CPointElemSurf<T> Nearest_Point_MeshTetFace3D(const delfem2::CVec3<T>& p0,
                                           const std::vector<double>& aXYZ,
                                           const std::vector<int>& aTet,
                                           const std::vector<int>& aTetFaceSrf);
  
template <typename T>
CPointElemSolid<T> Nearest_Point_MeshTet3D(const delfem2::CVec3<T>& q,
                                        const std::vector<double>& aXYZ,
                                        const std::vector<int>& aTet);
  
template <typename T>
CPointElemSolid<T> Nearest_Point_MeshTet3D(const delfem2::CVec3<T>& p,
                                        int itet_start, // starting triangle
                                        const std::vector<double>& aXYZ,
                                        const std::vector<int>& aTet,
                                        const std::vector<int>& aTetSurRel);

template <typename T>
double SDFNormal_NearestPoint(delfem2::CVec3<T>& n0,
                              const delfem2::CVec3<T>& p0,
                              const CPointElemSurf<T>& pes,
                              const std::vector<double>& aXYZ,
                              const std::vector<unsigned int>& aTri,
                              const std::vector<double>& aNorm);

template <typename T>
double SDFNormal_NearestPoint(delfem2::CVec3<T>& n0,
                              const delfem2::CVec3<T>& p0,
                              const CPointElemSurf<T>& pes,
                              const double* aXYZ, unsigned int nXYZ,
                              const unsigned int* aTri, unsigned int nTri,
                              const double* aNorm);

template <typename T>
double DistanceToTri(CPointElemSurf<T>& pes,
                     const delfem2::CVec3<T>& p,
                     unsigned int itri0,
                     const std::vector<double>& aXYZ,
                     const std::vector<unsigned int>& aTri);

template <typename T>
double DistanceToTri(CPointElemSurf<T>& pes,
                     const delfem2::CVec3<T>& p,
                     unsigned int itri0,
                     const double* aXYZ, unsigned int nXYZ,
                     const unsigned int* aTri, unsigned int nTri);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/srchuni_v3.cpp"
#endif
  
#endif /* search_mesh_hpp */
