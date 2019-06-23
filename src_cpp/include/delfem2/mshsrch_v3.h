#ifndef MSHSRCH_V3_H
#define MSHSRCH_V3_H

#include <stdio.h>

#include "delfem2/vec3.h"

////////////////////////////////////////////////////////////////////

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
  CVector3 getPos_Tet(const std::vector<double>& aXYZ, const std::vector<int>& aTet) const;
  void setPos_Tet(int it0, const CVector3& q, const std::vector<double>& aXYZ, const std::vector<int>& aTet);
public:
  int ielem;
  double r0, r1, r2;
};

// class for point in a triangle
class CPointElemSurf{
public:
  CPointElemSurf() : itri(-1), r0(0), r1(0) {}
  CPointElemSurf(int itri, double r0, double r1):itri(itri), r0(r0),r1(r1) {}
  CVector3 getPos_Tri(const std::vector<double>& aXYZ, const std::vector<int>& aTri) const;
  CVector3 getPos_TetFace(const std::vector<double>& aXYZ, const std::vector<int>& aTet, const std::vector<int>& aTetFace) const;
public:
  int itri;
  double r0, r1;
};
std::ostream &operator<<(std::ostream &output, const CPointElemSurf& v);
std::istream &operator>>(std::istream &input, CPointElemSurf& v);

////////////////////////////////////////////////////////////////////

CPointElemSurf intersect_Ray_Tri3D(double& depth,
                                   const CVector3& org, const CVector3& dir,
                                   int itri,
                                   const std::vector<unsigned int>& aTri,
                                   const std::vector<double>& aXYZ);
CPointElemSurf intersect_Ray_MeshTri3D(const CVector3& org, const CVector3& dir,
                                       const std::vector<unsigned int>& aTri,
                                       const std::vector<double>& aXYZ);
CPointElemSurf intersect_Ray_MeshTri3D(const CVector3& org, const CVector3& dir,
                                       int itri_start, // starting triangle
                                       const std::vector<unsigned int>& aTri,
                                       const std::vector<double>& aXYZ,
                                       const std::vector<int>& aTriSurRel);
CPointElemSurf intersect_Ray_MeshTri3D(const CVector3& org, const CVector3& dir,
                                       const CPointElemSurf& ptri0,
                                       const std::vector<unsigned int>& aTri,
                                       const std::vector<double>& aXYZ,
                                       const std::vector<int>& aTriSurRel);
CPointElemSurf intersect_Ray_MeshTriFlag3D(const CVector3& org, const CVector3& dir,
                                           const std::vector<unsigned int>& aTri,
                                           const std::vector<double>& aXYZ,
                                           int iflag,
                                           const std::vector<int>& aFlag);

/////////////////////////////////////////////////////////////////////////////////////////


CPointElemSurf nearest_Point_MeshTri3D(const CVector3& q,
                                       const std::vector<double>& aXYZ,
                                       const std::vector<int>& aTri);
CPointElemSurf nearest_Point_MeshTetFace3D(const CVector3& p0,
                                           const std::vector<double>& aXYZ,
                                           const std::vector<int>& aTet,
                                           const std::vector<int>& aTetFaceSrf);
CPointElemSolid nearest_Point_MeshTet3D(const CVector3& q,
                                        const std::vector<double>& aXYZ,
                                        const std::vector<int>& aTet);
CPointElemSolid nearest_Point_MeshTet3D(const CVector3& p,
                                        int itet_start, // starting triangle
                                        const std::vector<double>& aXYZ,
                                        const std::vector<int>& aTet,
                                        const std::vector<int>& aTetSurRel);

#endif /* search_mesh_hpp */
