/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_MSHSRCH_V3_H
#define DFM2_MSHSRCH_V3_H

#include <stdio.h>
#include <map>
#include "delfem2/vec3.h"

// -----------------------------------

namespace delfem2 {
  
bool intersectRay_Tri3D(double& r0, double& r1,
                        const CVec3& org, const CVec3& dir,
                        const CVec3& p0, const CVec3& p1, const CVec3& p2);

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
  delfem2::CVec3 getPos_Tet(const std::vector<double>& aXYZ, const std::vector<int>& aTet) const;
  void setPos_Tet(int it0, const delfem2::CVec3& q, const std::vector<double>& aXYZ, const std::vector<int>& aTet);
public:
  int ielem;
  double r0, r1, r2;
};

// class for point in a triangle
class CPointElemSurf{
public:
  CPointElemSurf() : itri(-1), r0(0), r1(0) {}
  CPointElemSurf(int itri, double r0, double r1):itri(itri), r0(r0),r1(r1) {}
  delfem2::CVec3 Pos_Tri(const std::vector<double>& aXYZ,
                   const std::vector<unsigned int>& aTri) const;
  delfem2::CVec3 Pos_Tri(const double* aXYZ, unsigned int nXYZ,
                   const unsigned int* aTri, unsigned int nTri) const;
  delfem2::CVec3 Pos_TetFace(const std::vector<double>& aXYZ,
                       const std::vector<int>& aTet,
                       const std::vector<int>& aTetFace) const;
  delfem2::CVec3 UNorm_Tri(const std::vector<double>& aXYZ,
                     const std::vector<unsigned int>& aTri,
                     const std::vector<double>& aNorm) const;
  delfem2::CVec3 UNorm_Tri(const double* aXYZ, unsigned int nXYZ,
                     const unsigned int* aTri, unsigned int nTri,
                     const double* aNorm) const;
  bool Check(const std::vector<double>& aXYZ,
             const std::vector<unsigned int>& aTri,
             double eps) const;
public:
  int itri; // can be -1
  double r0, r1;
};
std::ostream &operator<<(std::ostream &output, const CPointElemSurf& v);
std::istream &operator>>(std::istream &input, CPointElemSurf& v);

// ----------------------------------------------------------

std::vector<CPointElemSurf>
IntersectionLine_MeshTri3D(
    const delfem2::CVec3& org, const delfem2::CVec3& dir,
    const std::vector<unsigned int>& aTri,
    const std::vector<double>& aXYZ);

void IntersectionRay_MeshTri3D (
    std::map<double,CPointElemSurf>& mapDepthPES,
    const delfem2::CVec3& org, const delfem2::CVec3& dir,
    const std::vector<unsigned int>& aTri,
    const std::vector<double>& aXYZ);

void IntersectionRay_MeshTri3DPart(std::map<double,CPointElemSurf>& mapDepthPES,
                                   const delfem2::CVec3& org, const delfem2::CVec3& dir,
                                   const std::vector<unsigned int>& aTri,
                                   const std::vector<double>& aXYZ,
                                   const std::vector<int>& aIndTri);
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

// above functions for ray interesection
// -----------------------------------------------------------
// below functions for nearest


CPointElemSurf Nearest_Point_MeshTri3D(const delfem2::CVec3& q,
                                       const std::vector<double>& aXYZ,
                                       const std::vector<unsigned int>& aTri);
CPointElemSurf Nearest_Point_MeshTri3DPart(const delfem2::CVec3& q,
                                           const std::vector<double>& aXYZ,
                                           const std::vector<unsigned int>& aTri,
                                           const std::vector<int>& aIndTri_Cand);
CPointElemSurf Nearest_Point_MeshTetFace3D(const delfem2::CVec3& p0,
                                           const std::vector<double>& aXYZ,
                                           const std::vector<int>& aTet,
                                           const std::vector<int>& aTetFaceSrf);
CPointElemSolid Nearest_Point_MeshTet3D(const delfem2::CVec3& q,
                                        const std::vector<double>& aXYZ,
                                        const std::vector<int>& aTet);
CPointElemSolid Nearest_Point_MeshTet3D(const delfem2::CVec3& p,
                                        int itet_start, // starting triangle
                                        const std::vector<double>& aXYZ,
                                        const std::vector<int>& aTet,
                                        const std::vector<int>& aTetSurRel);

double SDFNormal_NearestPoint(delfem2::CVec3& n0,
                              const delfem2::CVec3& p0,
                              const CPointElemSurf& pes,
                              const std::vector<double>& aXYZ,
                              const std::vector<unsigned int>& aTri,
                              const std::vector<double>& aNorm);

double SDFNormal_NearestPoint(delfem2::CVec3& n0,
                              const delfem2::CVec3& p0,
                              const CPointElemSurf& pes,
                              const double* aXYZ, unsigned int nXYZ,
                              const unsigned int* aTri, unsigned int nTri,
                              const double* aNorm);

double DistanceToTri(CPointElemSurf& pes,
                     const delfem2::CVec3& p,
                     unsigned int itri0,
                     const std::vector<double>& aXYZ,
                     const std::vector<unsigned int>& aTri);

double DistanceToTri(CPointElemSurf& pes,
                     const delfem2::CVec3& p,
                     unsigned int itri0,
                     const double* aXYZ, unsigned int nXYZ,
                     const unsigned int* aTri, unsigned int nTri);

}
  
#endif /* search_mesh_hpp */
