#include <set>
#include <algorithm>

#include "delfem2/dyntri_v3.h"


CVector3 normalTri
(int itri0,
 const std::vector<CEPo2>& aPo3D,
 const std::vector<ETri>& aSTri,
 const std::vector<CVector3>& aXYZ)
{
  int i0 = aSTri[itri0].v[0];
  int i1 = aSTri[itri0].v[1];
  int i2 = aSTri[itri0].v[2];
  CVector3 n = Normal(aXYZ[i0], aXYZ[i1], aXYZ[i2]);
  return n.Normalize();
}


bool CheckTri
(const std::vector<CEPo2>& aPo3D,
 const std::vector<ETri>& aSTri,
 const std::vector<CVector3>& aXYZ,
 bool is_assert)
{
  for (int itri = 0; itri<aSTri.size(); itri++){
    const ETri& ref_tri = aSTri[itri];
    const int i0 = ref_tri.v[0];
    if( i0 == -1 ) continue;
    const int i1 = ref_tri.v[1];
    const int i2 = ref_tri.v[2];
    assert( i0 >=0 && i0 < (int)aPo3D.size() );
    assert( i1 >=0 && i1 < (int)aPo3D.size() );
    assert( i2 >=0 && i2 < (int)aPo3D.size() );
    double area = TriArea(aXYZ[i0], aXYZ[i1], aXYZ[i2]);
    if (area<1.0e-10){ // negative volume
    }
  }
  return true;
}

/*
template <typename TYPE>
void makeNormal
(std::vector<CEPo<TYPE>>& aPo3D,
const std::vector<STri2D>& aSTri)
{
  for (int ip = 0; ip<aPo3D.size(); ip++){
    aPo3D[ip].n.SetZero();
  }
  for (int itri = 0; itri<aSTri.size(); itri++){
    int i0 = aSTri[itri].v[0];
    int i1 = aSTri[itri].v[1];
    int i2 = aSTri[itri].v[2];
    CVector3 n = Normal(aPo3D[i0].p, aPo3D[i1].p, aPo3D[i2].p);
    aPo3D[i0].n += n;
    aPo3D[i1].n += n;
    aPo3D[i2].n += n;
  }
  for (int ip = 0; ip<aPo3D.size(); ip++){
    aPo3D[ip].n.SetNormalizedVector();
  }
}*/


/*
template <typename TYPE>
bool FindEdge
(unsigned int& itri0, unsigned int& inotri0, unsigned int& inotri1,
 ///
 const unsigned int& ipo0, const unsigned int& ipo1,  
 const std::vector<CEPo<TYPE>>& po, const std::vector<STri2D>& tri)
{
  const unsigned int itri_ini = po[ipo0].e;
  const unsigned int inotri_ini = po[ipo0].d;
  unsigned int inotri_cur = inotri_ini;
  unsigned int itri_cur = itri_ini;
  for (;;){	
    assert(tri[itri_cur].v[inotri_cur]==ipo0);
    {	
      const unsigned int inotri2 = indexRot3[1][inotri_cur];
      if (tri[itri_cur].v[inotri2]==ipo1){
        itri0 = itri_cur;
        inotri0 = inotri_cur;
        inotri1 = inotri2;
        assert(tri[itri0].v[inotri0]==ipo0);
        assert(tri[itri0].v[inotri1]==ipo1);
        return true;
      }
    }
    {	
      const unsigned int inotri2 = indexRot3[2][inotri_cur];
      if (tri[itri_cur].s2[inotri2]==-1){ break; }
      const unsigned int itri_nex = tri[itri_cur].s2[inotri2];
      const unsigned int* rel = relTriTri[tri[itri_cur].r2[inotri2]];
      const unsigned int inotri3 = rel[inotri_cur];
      assert(tri[itri_nex].v[inotri3]==ipo0);
      if (itri_nex==itri_ini) return false;
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
  }

  inotri_cur = inotri_ini;
  itri_cur = itri_ini;
  for (;;){	
    assert(tri[itri_cur].v[inotri_cur]==ipo0);
    {
      const unsigned int inotri2 = indexRot3[1][inotri_cur];
      if (tri[itri_cur].s2[inotri2]==-1){ break; }
      const unsigned int itri_nex = tri[itri_cur].s2[inotri2];
      const unsigned int* rel = relTriTri[tri[itri_cur].r2[inotri2]];
      const unsigned int inotri3 = rel[inotri_cur];
      assert(tri[itri_nex].v[inotri3]==ipo0);
      if (itri_nex==itri_ini){	// 
        itri0 = 0;
        inotri0 = 0; inotri1 = 0;
        return false;
      }
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
    {
      const unsigned int inotri2 = indexRot3[1][inotri_cur];
      if (tri[itri_cur].v[inotri2]==ipo1){
        itri0 = itri_cur;
        inotri0 = inotri_cur;
        inotri1 = inotri2;
        assert(tri[itri0].v[inotri0]==ipo0);
        assert(tri[itri0].v[inotri1]==ipo1);
        return true;
      }
    }
  }

  return false;
}
*/


int InsertPoint_Mesh
(const int itri0,
 double& r0,
 double& r1,
 std::vector<CEPo2>& aPo3D,
 std::vector<ETri>& aSTri,
 std::vector<CVector3>& aXYZ)
{
  if (itri0==-1) return -1;
  CVector3 pos,norm;
  {
    const int i0 = aSTri[itri0].v[0];
    const int i1 = aSTri[itri0].v[1];
    const int i2 = aSTri[itri0].v[2];
    const CVector3& p0 = aXYZ[i0];
    const CVector3& p1 = aXYZ[i1];
    const CVector3& p2 = aXYZ[i2];
    pos = r0*p0 + r1*p1 + (1-r0-r1)*p2;
    UnitNormal(norm, p0, p1, p2);
  }
  int ipo_ins = (int)aPo3D.size();
  aPo3D.push_back( CEPo2());
  aXYZ.push_back(pos);
  //  if( ptri.iedge == -1 ){ // inside tri
  InsertPoint_Elem(ipo_ins, itri0, aPo3D, aSTri);
  /*
   }
   else{
   InsertPoint_ElemEdge(ipo_ins, itri0, ptri.iedge, aPo3D, aSTri);
   }
   */
  return ipo_ins;
}


bool pickMesh
(CVector3& p,
 int& itriOut,
 double& r0Out,
 double& r1Out,
 ////
 const CVector3& org,
 const CVector3& dir,
 int itri_start, // starting triangle
 const std::vector<CEPo2>& aPo3D,
 const std::vector<ETri>& aSTri,
 std::vector<CVector3>& aXYZ)
{
  int itri1 = itri_start;
  for (int itr = 0; itr<50; itr++){
    if (itri1==-1) return false;
    int ip0 = aSTri[itri1].v[0];
    int ip1 = aSTri[itri1].v[1];
    int ip2 = aSTri[itri1].v[2];
    const CVector3& tp0 = aXYZ[ip0];
    const CVector3& p1 = aXYZ[ip1];
    const CVector3& p2 = aXYZ[ip2];
    const double v0 = volume_Tet(p1, p2, org, org+dir);
    const double v1 = volume_Tet(p2, tp0, org, org+dir);
    const double v2 = volume_Tet(tp0, p1, org, org+dir);
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


int pickTriangle
(CVector3& p,
 const CVector3& org, const CVector3& dir,
 int itri_start, // starting triangle
 const std::vector<CEPo2>& aPo3D,
 const std::vector<ETri>& aSTri,
 const std::vector<CVector3>& aXYZ)
{
  int itri1 = itri_start;
  for (int itr = 0; itr<50; itr++){
    if (itri1==-1) return -1;
    int ip0 = aSTri[itri1].v[0];
    int ip1 = aSTri[itri1].v[1];
    int ip2 = aSTri[itri1].v[2];
    const CVector3& tp0 = aXYZ[ip0];
    const CVector3& p1 = aXYZ[ip1];
    const CVector3& p2 = aXYZ[ip2];
    double v0 = volume_Tet(p1, p2, org, org+dir);
    double v1 = volume_Tet(p2, tp0, org, org+dir);
    double v2 = volume_Tet(tp0, p1, org, org+dir);
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


bool FindRayTriangleMeshIntersectionClosestToPoint
(const CVector3 &line0,
 const CVector3 &line1,
 const std::vector<ETri>& aTri,
 const std::vector<CEPo2> &aPoint3D,
 std::vector<CVector3>& aVec3,
 const CVector3 &targetPoint,
 CVector3 &intersectionPoint)
{
	std::vector<CVector3> intersectionPoints;
	if (!FindRayTriangleMeshIntersections(line0, line1,
                                        aTri, aPoint3D, aVec3,
                                        intersectionPoints))
	{
		return false;
	}
  
	// Find the point that is the closest to the target point
	float minSquareDistance = 1.0e16;
	for (unsigned int i = 0; i < intersectionPoints.size(); i++)
	{
		float currSquareDistance = 
    (intersectionPoints[i].x - targetPoint.x) * (intersectionPoints[i].x - targetPoint.x) +
    (intersectionPoints[i].y - targetPoint.y) * (intersectionPoints[i].y - targetPoint.y) +
    (intersectionPoints[i].z - targetPoint.z) * (intersectionPoints[i].z - targetPoint.z);
		if (currSquareDistance < minSquareDistance)
		{
			intersectionPoint = intersectionPoints[i];
			minSquareDistance = currSquareDistance;
		}
	}
  
	return true;
}


bool FindRayTriangleMeshIntersections
(const CVector3 &line0,
 const CVector3 &line1,
 const std::vector<ETri>& aTri,
 const std::vector<CEPo2> &aPoint3D,
 std::vector<CVector3>& aVec3,
 std::vector<CVector3> &intersectionPoints)
{
	intersectionPoints.clear();
  
	// Find all the intersection points between this ray and all triangles in the mesh
	for (unsigned int i = 0; i < aTri.size(); i++)
	{
		CVector3 intersectionPoint;
		if (isRayIntersectingTriangle(line0, line1,
                                  aVec3[aTri[i].v[0]],
                                  aVec3[aTri[i].v[1]],
                                  aVec3[aTri[i].v[2]],
                                  intersectionPoint))
		{
			intersectionPoints.push_back(intersectionPoint);
		}
	}  
	if (intersectionPoints.empty())
	{
		return false;
	} else {
		return true;
	}  
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

bool DelaunayAroundPoint
(int ipo0,
 std::vector<CEPo2>& aPo,
 std::vector<ETri>& aTri,
 const std::vector<CVector3>& aVec3)
{
  assert(ipo0 < (int)aPo.size());
  if (aPo[ipo0].e==-1) return true;
  
  assert(aPo[ipo0].e>=0&&(int)aPo[ipo0].e < (int)aTri.size());
  assert(aTri[aPo[ipo0].e].v[aPo[ipo0].d]==ipo0);
  
  const int itri0 = aPo[ipo0].e;
  int inotri0 = aPo[ipo0].d;
  
  int itri_cur = itri0;
  int inotri_cur = aPo[ipo0].d;
  bool flag_is_wall = false;
  for (;;){
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    if (aTri[itri_cur].s2[inotri_cur]>=0&&aTri[itri_cur].s2[inotri_cur]<(int)aTri.size()){
      assert(aTri[itri_cur].v[inotri_cur]==ipo0);
      // check opposing element
      const int itri_dia = aTri[itri_cur].s2[inotri_cur];
      const unsigned int* rel_dia = relTriTri[aTri[itri_cur].r2[inotri_cur]];
      const int inotri_dia = rel_dia[inotri_cur];
      assert(aTri[itri_dia].s2[inotri_dia]==itri_cur);
      const int ipo_dia = aTri[itri_dia].v[inotri_dia];
      if (DetDelaunay(aVec3[aTri[itri_cur].v[0]],
                      aVec3[aTri[itri_cur].v[1]],
                      aVec3[aTri[itri_cur].v[2]],
                      aVec3[ipo_dia])==0)
      {
        bool res = FlipEdge(itri_cur, inotri_cur, aPo, aTri);
        if( res ){
          inotri_cur = 2;
          assert(aTri[itri_cur].v[inotri_cur]==ipo0);
          if (itri_cur==itri0) inotri0 = inotri_cur;
          continue;
        }
        else{
          break;
        }
      }
    }
    MoveCCW(itri_cur, inotri_cur, flag_is_wall, aPo,aTri);
    if( flag_is_wall ) break;
    if( itri_cur == itri0 ) break;
    /*
     {  // next element
     const int inotri1 = indexRot3[1][inotri_cur];
     if (aTri[itri_cur].s2[inotri1]==-1){
     flag_is_wall = true;
     break;
     }
     const int itri_nex = aTri[itri_cur].s2[inotri1];
     const unsigned int* rel_nex = relTriTri[aTri[itri_cur].r2[inotri1]];
     const int inotri_nex = rel_nex[inotri_cur];
     assert(aTri[itri_nex].v[inotri_nex]==ipo0);
     if (itri_nex==itri0) break;  // finish if we reach starting elemnt
     itri_cur = itri_nex;
     inotri_cur = inotri_nex;
     }
     */
  }
  if (!flag_is_wall) return true;
  
  ////////////////////////////////
  // rotate counter clock-wise
  
  itri_cur = itri0;
  inotri_cur = inotri0;
  for (;;){
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    
    if (aTri[itri_cur].s2[inotri_cur]>=0&&aTri[itri_cur].s2[inotri_cur]<(int)aTri.size()){
      // check elements in opposing side
      const int itri_dia = aTri[itri_cur].s2[inotri_cur];
      const unsigned int* rel_dia = relTriTri[aTri[itri_cur].r2[inotri_cur]];
      const int inotri_dia = rel_dia[inotri_cur];
      assert(aTri[itri_dia].s2[inotri_dia]==itri_cur);
      const int ipo_dia = aTri[itri_dia].v[inotri_dia];
      if (DetDelaunay(aVec3[aTri[itri_cur].v[0]],
                      aVec3[aTri[itri_cur].v[1]],
                      aVec3[aTri[itri_cur].v[2]],
                      aVec3[ipo_dia])==0)  // Delaunay condition is not satisfiled
      {
        FlipEdge(itri_cur, inotri_cur, aPo, aTri);
        itri_cur = itri_dia;
        inotri_cur = 1;
        assert(aTri[itri_cur].v[inotri_cur]==ipo0);
        continue;
      }
    }
    
    {
      const int inotri2 = (inotri_cur+2)%3; //  indexRot3[2][inotri_cur];
      if (aTri[itri_cur].s2[inotri2]==-1){
        return true;
      }
      const int itri_nex = aTri[itri_cur].s2[inotri2];
      const unsigned int* rel_nex = relTriTri[aTri[itri_cur].r2[inotri2]];
      const int inotri_nex = rel_nex[inotri_cur];
      assert(aTri[itri_nex].v[inotri_nex]==ipo0);
      assert(itri_nex!=itri0);  // finsih if reach starting elemnet
      itri_cur = itri_nex;
      inotri_cur = inotri_nex;
    }
  }
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
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






#ifdef USE_GL

#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/gl.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#else
#include <GL/gl.h>
#endif

void DrawMeshDynTri_FaceNorm
(const std::vector< CEPo2>& aPo3D,
 const std::vector<ETri>& aSTri,
 const std::vector<CVector3>& aVec3)
{
  //  ::glPushAttrib(GL_ENABLE_BIT);
  
  ::glEnable(GL_LIGHTING);
  //  ::myGlColorDiffuse(CColor::Orange());
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri=0; itri<aSTri.size(); ++itri){
    const int i0 = aSTri[itri].v[0];
    const int i1 = aSTri[itri].v[1];
    const int i2 = aSTri[itri].v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
      continue;
    }
    {
      CVector3 n; UnitNormal(n, aVec3[i0], aVec3[i1], aVec3[i2]);
      ::glNormal3d(n.x,n.y,n.z);
    }
    {
      CVector3 p0 = aVec3[i0];
      ::glVertex3d(p0.x,p0.y,p0.z);
    }
    {
      CVector3 p1 = aVec3[i1];
      ::glVertex3d(p1.x,p1.y,p1.z);
    }
    {
      CVector3 p2 = aVec3[i2];
      ::glVertex3d(p2.x,p2.y,p2.z);
    }
  }
  ::glEnd();
}

void DrawMeshDynTri_Edge
(const std::vector< CEPo2>& aPo3D,
 const std::vector<ETri>& aSTri,
 const std::vector<CVector3>& aVec3)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (unsigned int itri = 0; itri<aSTri.size(); ++itri){
    const int i0 = aSTri[itri].v[0];
    const int i1 = aSTri[itri].v[1];
    const int i2 = aSTri[itri].v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
    }
    const CVector3& p0 = aVec3[i0];
    const CVector3& p1 = aVec3[i1];
    const CVector3& p2 = aVec3[i2];
    glVertex3d(p0.x,p0.y,p0.z);
    glVertex3d(p1.x,p1.y,p1.z);
    
    glVertex3d(p1.x,p1.y,p1.z);
    glVertex3d(p2.x,p2.y,p2.z);
    
    glVertex3d(p2.x,p2.y,p2.z);
    glVertex3d(p0.x,p0.y,p0.z);
  }
  ::glEnd();
}

#endif
