/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <algorithm>
#include <climits>
#include "delfem2/geoproximity3_v3.h"
#include "delfem2/geodelaunay3_v3.h"
#include "delfem2/dtri3_v3dtri.h"

//namespace dfm2 = delfem2;

// ----------------------------------------

namespace delfem2{
namespace dtriv3{

//! Volume of a tetrahedra
template <typename T>
T Volume_Tet(
    const CVec3<T>& v0,
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3 )
{
//  return delfem2::Volume_Tet3(v0.p, v1.p, v2.p, v3.p);
  double v = (v1.p[0]-v0.p[0])*( (v2.p[1]-v0.p[1])*(v3.p[2]-v0.p[2]) - (v3.p[1]-v0.p[1])*(v2.p[2]-v0.p[2]) )
             + (v1.p[1]-v0.p[1])*( (v2.p[2]-v0.p[2])*(v3.p[0]-v0.p[0]) - (v3.p[2]-v0.p[2])*(v2.p[0]-v0.p[0]) )
             + (v1.p[2]-v0.p[2])*( (v2.p[0]-v0.p[0])*(v3.p[1]-v0.p[1]) - (v3.p[0]-v0.p[0])*(v2.p[1]-v0.p[1]) );
  return v*0.16666666666666666666666666666667;
}

DFM2_INLINE int InsertPoint_Mesh(
    const int itri0,
    double& r0,
    double& r1,
    std::vector<CDynPntSur>& aPo3D,
    std::vector<CDynTri>& aSTri,
    std::vector<CVec3d>& aXYZ)
{
  if (itri0==-1) return -1;
  CVec3d pos;
  {
    const int i0 = aSTri[itri0].v[0];
    const int i1 = aSTri[itri0].v[1];
    const int i2 = aSTri[itri0].v[2];
    const CVec3d& p0 = aXYZ[i0];
    const CVec3d& p1 = aXYZ[i1];
    const CVec3d& p2 = aXYZ[i2];
    pos = r0*p0 + r1*p1 + (1-r0-r1)*p2;
  }
  int ipo_ins = (int)aPo3D.size();
  aPo3D.emplace_back( );
  aXYZ.push_back(pos);
  InsertPoint_Elem(ipo_ins, itri0, aPo3D, aSTri);
  return ipo_ins;
}


DFM2_INLINE bool pickMesh
 (CVec3d& p,
  int& itriOut,
  double& r0Out,
  double& r1Out,
  //
  const CVec3d& org,
  const CVec3d& dir,
  int itri_start, // starting triangle
  const std::vector<CDynTri>& aSTri,
  const std::vector<CVec3d>& aXYZ)
{
  int itri1 = itri_start;
  for (int itr = 0; itr<50; itr++){
    if (itri1==-1) return false;
    int ip0 = aSTri[itri1].v[0];
    int ip1 = aSTri[itri1].v[1];
    int ip2 = aSTri[itri1].v[2];
    const CVec3d& tp0 = aXYZ[ip0];
    const CVec3d& p1 = aXYZ[ip1];
    const CVec3d& p2 = aXYZ[ip2];
    const double v0 = Volume_Tet(p1, p2, org, org+dir);
    const double v1 = Volume_Tet(p2, tp0, org, org+dir);
    const double v2 = Volume_Tet(tp0, p1, org, org+dir);
    const double r0 = v0/(v0+v1+v2);
    const double r1 = v1/(v0+v1+v2);
    const double r2 = v2/(v0+v1+v2);
    /*
     if(      fabs(r0) < 1.0e-3 && r1>0 && r2>0 ){
     p = tp0*r0+p1*r1+p2*r2;
     CPointTri ptri;
     ptri.itri = itri1;
     ptri.iedge = 0;
     ptri.r0 = v1/(v1+v2);
     return true;
     }
     else if( fabs(r1) < 1.0e-3 && r2>0 && r0>0 ){
     p = tp0*r0+p1*r1+p2*r2;
     CPointTri ptri;
     ptri.itri = itri1;
     ptri.iedge = 1;
     ptri.r0 = v2/(v2+v0);
     return true;
     }
     else if( fabs(r2) < 1.0e-3 && r0>0 && r1>0 ){
     p = tp0*r0+p1*r1+p2*r2;
     CPointTri ptri;
     ptri.itri = itri1;
     ptri.iedge = 2;
     ptri.r0 = v0/(v0+v1);
     return true;
     }
     */
      ////
    if (r0>0&&r1>0&&r2>0){
      p = tp0*r0+p1*r1+p2*r2;
      itriOut = itri1;
      r0Out = r0;
      r1Out = r1;
      return true;
    }
    if (     r0<r1 && r0<r2){ itri1 = aSTri[itri1].s2[0]; }
    else if (r1<r0 && r1<r2){ itri1 = aSTri[itri1].s2[1]; }
    else{                     itri1 = aSTri[itri1].s2[2]; }
  }
  return false;
}


DFM2_INLINE int pickTriangle
 (CVec3d& p,
  const CVec3d& org,
  const CVec3d& dir,
  int itri_start, // starting triangle
  const std::vector<CDynTri>& aSTri,
  const std::vector<CVec3d>& aXYZ)
{
  int itri1 = itri_start;
  for (int itr = 0; itr<50; itr++){
    if (itri1==-1) return -1;
    int ip0 = aSTri[itri1].v[0];
    int ip1 = aSTri[itri1].v[1];
    int ip2 = aSTri[itri1].v[2];
    const CVec3d& tp0 = aXYZ[ip0];
    const CVec3d& p1 = aXYZ[ip1];
    const CVec3d& p2 = aXYZ[ip2];
    double v0 = Volume_Tet(p1, p2, org, org+dir);
    double v1 = Volume_Tet(p2, tp0, org, org+dir);
    double v2 = Volume_Tet(tp0, p1, org, org+dir);
    if (v0>0&&v1>0&&v2>0){
      double r0 = v0/(v0+v1+v2);
      double r1 = v1/(v0+v1+v2);
      double r2 = v2/(v0+v1+v2);
      p = tp0*r0+p1*r1+p2*r2;
      return itri1;
    }
    if (v0<v1 && v0<v2){
      itri1 = aSTri[itri1].s2[0];
    }
    else if (v1<v0 && v1<v2){
      itri1 = aSTri[itri1].s2[1];
    }
    else{
      itri1 = aSTri[itri1].s2[2];
    }
  }
  return -1;
}


DFM2_INLINE bool FindRayTriangleMeshIntersectionClosestToPoint
 (CVec3d &intersectionPoint,
  const CVec3d &line0,
  const CVec3d &line1,
  const std::vector<CDynTri>& aTri,
  const std::vector<CVec3d>& aVec3,
  const CVec3d &targetPoint)
{
  std::vector<CVec3d> intersectionPoints;
  if (!FindRayTriangleMeshIntersections(intersectionPoints,
                                        line0, line1,
                                        aTri, aVec3))
  {
    return false;
  }
  
    // Find the point that is the closest to the target point
  double minSquareDistance = 1.0e16;
  for (auto & i : intersectionPoints)
  {
    float currSquareDistance =
    (i.x - targetPoint.x) * (i.x - targetPoint.x) +
    (i.y - targetPoint.y) * (i.y - targetPoint.y) +
    (i.z - targetPoint.z) * (i.z - targetPoint.z);
    if (currSquareDistance < minSquareDistance)
    {
      intersectionPoint = i;
      minSquareDistance = currSquareDistance;
    }
  }
  
  return true;
}

} // namespace dtriv3
} // namespace delfem2

// static functions
// --------------------------------------------------------------------------
// expose functions

DFM2_INLINE delfem2::CVec3d delfem2::UnitNormal_DTri3(
    unsigned int itri0,
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec3d>& aP3)
{
  const unsigned int i0 = aSTri[itri0].v[0];
  const unsigned int i1 = aSTri[itri0].v[1];
  const unsigned int i2 = aSTri[itri0].v[2];
  const CVec3d n = Normal(aP3[i0], aP3[i1], aP3[i2]);
  return n.normalized();
}


bool delfem2::AssertMeshDTri2(
	[[maybe_unused]] const std::vector<CDynPntSur>& aPo3D,
	const std::vector<CDynTri>& aSTri,
	const std::vector<CVec3d>& aXYZ)
{
  for (const auto & ref_tri : aSTri){
    const int i0 = ref_tri.v[0];
    if( i0 == -1 ) continue;
    const int i1 = ref_tri.v[1];
    const int i2 = ref_tri.v[2];
    assert( i0 >=0 && i0 < (int)aPo3D.size() );
    assert( i1 >=0 && i1 < (int)aPo3D.size() );
    assert( i2 >=0 && i2 < (int)aPo3D.size() );
    double area = Area_Tri(aXYZ[i0], aXYZ[i1], aXYZ[i2]);
    if (area<1.0e-10){ // very small volume
      assert(0);
      abort();
    }
  }
  return true;
}



DFM2_INLINE bool delfem2::FindRayTriangleMeshIntersections(
    std::vector<CVec3d> &intersectionPoints,
    const CVec3d &line0,
    const CVec3d &line1,
    const std::vector<CDynTri>& aTri,
    const std::vector<CVec3d>& aVec3)
{
	intersectionPoints.clear();
  
	// Find all the intersection points between this ray and all triangles in the mesh
	for (const auto & i : aTri)
	{
		CVec3d intersectionPoint;
		if (isRayIntersectingTriangle(line0, line1,
                                  aVec3[i.v[0]],
                                  aVec3[i.v[1]],
                                  aVec3[i.v[2]],
                                  intersectionPoint))
		{
			intersectionPoints.push_back(intersectionPoint);
		}
	}
  return !intersectionPoints.empty();
}

// -----------------------------------------------------------------

DFM2_INLINE bool delfem2::DelaunayAroundPoint(
    const unsigned int ipo0,
    std::vector<CDynPntSur>& aPo,
    std::vector<CDynTri>& aTri,
    const std::vector<CVec3d>& aVec3)
{
  assert( ipo0 < aPo.size() );
  if ( aPo[ipo0].e==UINT_MAX ) return true;
  
  assert( aPo[ipo0].e<aTri.size() );
  assert( aTri[aPo[ipo0].e].v[aPo[ipo0].d]==ipo0 );
  
  unsigned int it0 = aPo[ipo0].e;
  unsigned int in0 = aPo[ipo0].d;
  bool flag_is_wall = false;
  for (;;){
    assert(aTri[it0].v[in0]==ipo0);
    if ( aTri[it0].s2[in0]<aTri.size() ){
      assert(aTri[it0].v[in0]==ipo0);
      const unsigned int jt0 = aTri[it0].s2[in0];
      const unsigned int jn0 = FindAdjEdgeIndex(aTri[it0],in0,aTri);
      assert(aTri[jt0].s2[jn0]==it0);
      const unsigned int jp0 = aTri[jt0].v[jn0];
      const int ires = DetDelaunay(
          aVec3[aTri[it0].v[0]],
          aVec3[aTri[it0].v[1]],
          aVec3[aTri[it0].v[2]],
          aVec3[jp0]);
      if( ires == 0 ){
        FlipEdge(it0, in0, aPo, aTri);
        in0 = 2;
        assert(aTri[it0].v[in0]==ipo0);
        continue;
      }
    }
    if( !MoveCCW(it0, in0, UINT_MAX, aTri) ){ flag_is_wall = true; break; }
    if( it0 == aPo[ipo0].e ) break;
  }
  if (!flag_is_wall) return true;
  
  // move couner-clock-wise
  // ----------------------------------
  // move clock-wise
  
  it0 = aPo[ipo0].e;
  in0 = aPo[ipo0].d;
  for (;;){
    assert( aTri[it0].v[in0]==ipo0 );
    if ( aTri[it0].s2[in0]<aTri.size() ){
      const unsigned int jt0 = aTri[it0].s2[in0];
      const unsigned int jn0 = FindAdjEdgeIndex(aTri[it0],in0,aTri);
      assert( aTri[jt0].s2[jn0]==it0 );
      const unsigned int jp0 = aTri[jt0].v[jn0];
      int ires = DetDelaunay(aVec3[aTri[it0].v[0]],
                             aVec3[aTri[it0].v[1]],
                             aVec3[aTri[it0].v[2]],
                             aVec3[jp0]);
      if( ires == 0 ){ // Delaunay condition is not satisfiled
        FlipEdge(it0, in0, aPo, aTri);
        it0 = jt0;
        in0 = 1;
        assert(aTri[it0].v[in0]==ipo0);
        continue;
      }
    }
    if( !MoveCW(it0, in0, UINT_MAX, aTri) ){ return  true; }
  }
  return true;
}

// ------------------------------------
/*
template <typename TYPE>
void GetTriAryAroundPoint
(unsigned int ipoin,
 std::vector< std::pair<unsigned int,unsigned int> >& aTriSurPo,
 const std::vector<CEPo<TYPE>>& aPo, const std::vector<STri2D>& aTri)
{
	const unsigned int itri_ini = aPo[ipoin].e;
	const unsigned int inoel_c_ini = aPo[ipoin].d;
	assert( itri_ini < aTri.size() );
	assert( inoel_c_ini < 3 );
	assert( aTri[itri_ini].v[inoel_c_ini] == ipoin );
	unsigned int itri0= itri_ini;
	unsigned int inoel_c0 = inoel_c_ini;
	unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
	for(;;){
		assert( itri0 < aTri.size() );
		assert( inoel_c0 < 3 );
		assert( aTri[itri0].v[inoel_c0] == ipoin );
    aTriSurPo.push_back( std::make_pair(itri0,inoel_c0) );
    /////
//    std::cout << ipoin << " " << itri0 << " " << inoel_b0 << " " << aTri.size() << " " << aTri[itri0].s2[inoel_b0] << std::endl;
    if( aTri[itri0].s2[inoel_b0] == -1 ){ break; }
		assert( aTri[itri0].s2[inoel_b0] >=0 && aTri[itri0].s2[inoel_b0] < (int)aTri.size() );
		unsigned int itri1 = aTri[itri0].s2[inoel_b0];
		const unsigned int rel01 = aTri[itri0].r2[inoel_b0];
		const unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
		const unsigned int inoel_b1 = relTriTri[rel01][ noelTriEdge[inoel_c0][1] ];
		assert( itri1 < aTri.size() );
		assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );
		assert( aTri[itri1].v[inoel_c1] == ipoin );
		if( itri1 == itri_ini ) return;
		itri0 = itri1;
		inoel_c0 = inoel_c1;
		inoel_b0 = inoel_b1;
	}
}*/




