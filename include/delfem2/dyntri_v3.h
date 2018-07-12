#ifndef SURFACE_MESH_H
#define SURFACE_MESH_H

#include <map>
#include <algorithm>
#include <stack>

#include "delfem2/vec3.h"

//////////////////////////////////////////////////////////////////////////////////////////////////

static const unsigned int relTriTri[3][3] = {
	{ 0, 2, 1 }, //  0
	{ 2, 1, 0 }, //  1 
	{ 1, 0, 2 }, //  2
};

static const unsigned int invRelTriTri[3] = { 0, 1, 2 };

class ETri{
public:
	int v[3];	//
	int s2[3];	//!< index of face that is adjacent to ith edge; The 0th edge is the edge facing 0th vertex
  int r2[3];	//!< relationship of vertex index between two adjacent faces
};

template <typename TYPE>
class CEPo{
public:
	CEPo(){
    e = -1;
    d = 0;
  }
	CEPo( const CEPo& rhs )
		: e(rhs.e), d(rhs.d), p(rhs.p), n(rhs.n), t(rhs.t){}
  CEPo(const CVector3 p, int ielem, unsigned int idir)
    : p(p), e(ielem), d(idir){}
  CEPo(const CVector3& p, const CVector3& n) :
    p(p), n(n), e(-1), d(0){}
public:
  int e;  
  int d; 
  ////
  CVector3 p; // position
  CVector3 n; // normal
  /////
  TYPE t;
};

bool MakePointSurTri
( const std::vector<ETri>& aTri, const unsigned int npoin,
 unsigned int* const elsup_ind, unsigned int& nelsup, unsigned int*& elsup );

void MakeEdge
(unsigned int* const edge_ind, unsigned int& nedge, unsigned int*& edge,
const std::vector<ETri>& aTri, const unsigned int npoin,
const unsigned int* elsup_ind, const unsigned int nelsup, const unsigned int* elsup);

bool MakeInnerRelationTri
(std::vector<ETri>& aTri, const unsigned int npoin,
const unsigned int* elsup_ind, const unsigned int nelsup, const unsigned int* elsup);

CVector3 ProjectPointOnTriangle
(const CVector3 &p0,
const CVector3 &tri_p1, const CVector3 &tri_p2, const CVector3 &tri_p3);

bool isPointInsideTriangle(const CVector3 &p0,
  const CVector3 &tri_p1, const CVector3 &tri_p2, const CVector3 &tri_p3);
bool isPointSameSide(const CVector3 &p0, const CVector3 &p1,
  const CVector3 &line_p0, const CVector3 &line_p1);
bool isRayIntersectingTriangle(const CVector3 &line0, const CVector3 &line1,
  const CVector3 &tri0, const CVector3 &tri1, const CVector3 &tri2,
  CVector3 &intersectionPoint);

bool CheckTri(const std::vector<ETri>& aTri);

template <typename TYPE>
bool CheckTri
(const std::vector< CEPo<TYPE> >& aPo3D,
 const std::vector<ETri>& aSTri,
 bool is_assert=true)
{
  const int npo = (int)aPo3D.size();
  const int ntri = (int)aSTri.size();

  for (int itri = 0; itri<ntri; itri++){
    const ETri& tri0 = aSTri[itri];
    if( tri0.v[0] == -1 ){
      assert(tri0.v[1] == -1);
      assert(tri0.v[2] == -1);
      continue;
    }
    if (tri0.v[0]==tri0.v[1]){ //     assert(tri0.v[2]!=tri0.v[0]);
      if( is_assert ){ assert(0); }
      return false;
    }
    if (tri0.v[1]==tri0.v[2]){ //     assert(tri0.v[2]!=tri0.v[0]);
      if( is_assert ){ assert(0); }
      return false;
    }
    if (tri0.v[2]==tri0.v[0]){ //     assert(tri0.v[2]!=tri0.v[0]);
      if( is_assert ){ assert(0); }
      return false;
    }
    ////
    if( tri0.s2[0]==tri0.s2[1] ){ // assert(tri0.s2[0]!=tri0.s2[2]);
//      std::cout << tri0.s2[0] << " " << tri0.s2[1] << std::endl;
      if( is_assert ){ assert(0); }
      return false;
    }
    if( tri0.s2[1]==tri0.s2[2] ){ // assert(tri0.s2[1]!=tri0.s2[2]);
      if( is_assert ){ assert(0); }
      return false;
    }
    if( tri0.s2[2]==tri0.s2[0] ){ // assert(tri0.s2[2]!=tri0.s2[0]);
      if( is_assert ){ assert(0); }
      return false;
    }
    /////
    for (unsigned int inotri = 0; inotri<3; inotri++){
      assert(tri0.v[inotri] < npo);
    }
    for (int iedtri = 0; iedtri<3; iedtri++){
      if (tri0.s2[iedtri]>=0&&tri0.s2[iedtri]<ntri){
        const int itri_s = tri0.s2[iedtri];
        const int irel = tri0.r2[iedtri];
        assert(itri_s < ntri);
        assert(irel < 3);
        {
          const int noel_dia = relTriTri[irel][iedtri];
          assert(noel_dia < 3);
          if (aSTri[itri_s].s2[noel_dia]!=itri){ // neibough relation broken
//            std::cout << itri << " " << itri_s << " " << noel_dia << std::endl;
          }
          assert(aSTri[itri_s].s2[noel_dia]==itri);
        }
        // check relation 
        for (int inoed = 0; inoed<2; inoed++){
          const int inoel = (iedtri+1+inoed)%3;//noelTriEdge[iedtri][inoed];
          if (tri0.v[inoel]!=aSTri[itri_s].v[(int)relTriTri[irel][inoel]]){
          }
          assert(tri0.v[inoel]==aSTri[itri_s].v[(int)relTriTri[irel][inoel]]);
        }
      }
    }
  }
  for (int ipoin = 0; ipoin<npo; ++ipoin){
    const int itri0 = aPo3D[ipoin].e;
    const int inoel0 = aPo3D[ipoin].d;
    if (aPo3D[ipoin].e>=0){
      assert(aPo3D[ipoin].d>=0&&aPo3D[ipoin].d < 3);
      if (aSTri[itri0].v[inoel0]!=ipoin){}
      assert(aSTri[itri0].v[inoel0]==ipoin);
    }
  }

  for (int itri = 0; itri<ntri; itri++){
    const ETri& ref_tri = aSTri[itri];
    const int i0 = ref_tri.v[0];
    if( i0 == -1 ) continue;
    const int i1 = ref_tri.v[1];
    const int i2 = ref_tri.v[2];
    assert( i0 >=0 && i0 < (int)aPo3D.size() );
    assert( i1 >=0 && i1 < (int)aPo3D.size() );
    assert( i2 >=0 && i2 < (int)aPo3D.size() );
    double area = TriArea(aPo3D[i0].p, aPo3D[i1].p, aPo3D[i2].p);
    if (area<1.0e-10){ // negative volume
    }
  }
  return true;
}


template <typename TYPE>
bool FindEdge
(int& itri0, int& inotri0, int& inotri1,
///
const int ipo0, const int ipo1,
const std::vector< CEPo<TYPE> >& aPo, const std::vector<ETri>& aTri)
{
  const int itri_ini = aPo[ipo0].e;
  const int inotri_ini = aPo[ipo0].d;
  int inotri_cur = inotri_ini;
  int itri_cur = itri_ini;
  for (;;){	// serch clock-wise
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    {
      const int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      if (aTri[itri_cur].v[inotri2]==ipo1){
        itri0 = itri_cur;
        inotri0 = inotri_cur;
        inotri1 = inotri2;
        assert(aTri[itri0].v[inotri0]==ipo0);
        assert(aTri[itri0].v[inotri1]==ipo1);
        return true;
      }
    }
    {
      const int inotri2 = (inotri_cur+2)%3; // indexRot3[2][inotri_cur];
      if (aTri[itri_cur].s2[inotri2]==-1){ break; }
      const int itri_nex = aTri[itri_cur].s2[inotri2];
      const unsigned int* rel = relTriTri[aTri[itri_cur].r2[inotri2]];
      const int inotri3 = rel[inotri_cur];
      assert(aTri[itri_nex].v[inotri3]==ipo0);
      if (itri_nex==itri_ini) return false;
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
  }

  inotri_cur = inotri_ini;
  itri_cur = itri_ini;
  for (;;){
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    {
      const int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      if (aTri[itri_cur].s2[inotri2]==-1){ break; }
      const int itri_nex = aTri[itri_cur].s2[inotri2];
      const unsigned int* rel = relTriTri[aTri[itri_cur].r2[inotri2]];
      const int inotri3 = rel[inotri_cur];
      assert(aTri[itri_nex].v[inotri3]==ipo0);
      if (itri_nex==itri_ini){	// end if it goes around
        itri0 = 0;
        inotri0 = 0; inotri1 = 0;
        return false;
      }
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
    {
      const int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      if (aTri[itri_cur].v[inotri2]==ipo1){
        itri0 = itri_cur;
        inotri0 = inotri_cur;
        inotri1 = inotri2;
        assert(aTri[itri0].v[inotri0]==ipo0);
        assert(aTri[itri0].v[inotri1]==ipo1);
        return true;
      }
    }
  }
  return false;
}

template <typename TYPE>
void makeNormal(
    std::vector<CEPo<TYPE> >& aPo3D,
    const std::vector<ETri>& aSTri)
{
  for (int ip = 0; ip<(int)aPo3D.size(); ip++){
    aPo3D[ip].n.SetZero();
  }
  for (int itri = 0; itri<(int)aSTri.size(); itri++){
    const int i0 = aSTri[itri].v[0];
    const int i1 = aSTri[itri].v[1];
    const int i2 = aSTri[itri].v[2];
    if( i0 == -1 ){
      assert(i1==-1); assert(i2==-1);
      continue;
    }
    CVector3 n = Normal(aPo3D[i0].p, aPo3D[i1].p, aPo3D[i2].p);
    aPo3D[i0].n += n;
    aPo3D[i1].n += n;
    aPo3D[i2].n += n;
  }
  for (int ip = 0; ip<(int)aPo3D.size(); ip++){
    if( aPo3D[ip].e == -1 ){
      aPo3D[ip].n = CVector3(1,0,0);
      continue;
    }
    aPo3D[ip].n.SetNormalizedVector();
  }
}

template <typename TYPE>
CVector3 normalTri
(int itri0,
 const std::vector<CEPo<TYPE> >& aPo3D,
 const std::vector<ETri>& aSTri)
{
  int i0 = aSTri[itri0].v[0];
  int i1 = aSTri[itri0].v[1];
  int i2 = aSTri[itri0].v[2];
  CVector3 n = Normal(aPo3D[i0].p, aPo3D[i1].p, aPo3D[i2].p);
  return n.Normalize();
}

///////////////////////////////////////////////////////////////
// topology edit

template <typename TYPE>
bool InsertPoint_ElemEdge
(const int ipo_ins,    //the index of the new point
 const int itri_ins,  //triangle index				
 const int ied_ins,  //edge index
 std::vector<CEPo<TYPE> >& aPo,
 std::vector<ETri>& aTri )
{
  assert(itri_ins < (int)aTri.size());
  assert(ipo_ins < (int)aPo.size());

  if (aTri[itri_ins].s2[ied_ins]==-1){
    assert(0);
  }
  
  // (node index opp to 0)*3+(node index opp to 1) -> relation index
  const int noel2RelTriTri[9] = {
    -1,	// 0 00
    -1,	// 1 01
    0,	// 2 02
    2, // 3 10
    -1, // 4 11
    -1,	// 5 12
    -1,	// 6 20
    1,	// 7 21
    -1, // 8 22
  };

  const int itri_adj = aTri[itri_ins].s2[ied_ins];
  const int ied_adj = (int)relTriTri[(int)aTri[itri_ins].r2[ied_ins]][ied_ins];
  assert(itri_adj < (int)aTri.size());
  assert(ied_ins < 3);

  const int itri0 = itri_ins;
  const int itri1 = itri_adj;
  const int itri2 = (int)aTri.size();
  const int itri3 = (int)aTri.size()+1;

  aTri.resize(aTri.size()+2);

  ETri old0 = aTri[itri_ins];
  ETri old1 = aTri[itri_adj];

  const int ino0_0 = ied_ins;
  const int ino1_0 = (ied_ins+1)%3; //noelTriEdge[ied_ins][0];
  const int ino2_0 = (ied_ins+2)%3; //noelTriEdge[ied_ins][1];

  const int ino0_1 = ied_adj;
  const int ino1_1 = (ied_adj+1)%3; //noelTriEdge[ied_adj][0];
  const int ino2_1 = (ied_adj+2)%3; //noelTriEdge[ied_adj][1];

  assert(old0.v[ino1_0]==old1.v[ino2_1]);
  assert(old0.v[ino2_0]==old1.v[ino1_1]);
  assert(old0.s2[ino0_0]==itri1);
  assert(old1.s2[ino0_1]==itri0);

  aPo[ipo_ins].e = itri0;			aPo[ipo_ins].d = 0;
  aPo[old0.v[ino2_0]].e = itri0;	aPo[old0.v[ino2_0]].d = 1;
  aPo[old0.v[ino0_0]].e = itri1;	aPo[old0.v[ino0_0]].d = 1;
  aPo[old1.v[ino2_1]].e = itri2;	aPo[old1.v[ino2_1]].d = 1;
  aPo[old1.v[ino0_1]].e = itri3;	aPo[old1.v[ino0_1]].d = 1;

  {
    ETri& ref_tri = aTri[itri0];
    ////////////////
    ref_tri.v[0] = ipo_ins;			    ref_tri.v[1] = old0.v[ino2_0];	ref_tri.v[2] = old0.v[ino0_0];
    ref_tri.s2[0] = old0.s2[ino1_0];	ref_tri.s2[1] = itri1;			    ref_tri.s2[2] = itri3;
    ////////////////
    if (old1.s2[ino1_0]>=0&&old1.s2[ino1_0]<(int)aTri.size()){
      assert(old0.r2[ino1_0] < 3);
      const unsigned int* rel = relTriTri[old0.r2[ino1_0]];
      ref_tri.r2[0] = noel2RelTriTri[rel[ino1_0]*3+rel[ino2_0]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old0.s2[ino1_0] < (int)aTri.size());
      aTri[old0.s2[ino1_0]].s2[rel[ino1_0]] = itri0;
      aTri[old0.s2[ino1_0]].r2[rel[ino1_0]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri1];
    ////////////////
    ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old0.v[ino0_0];	ref_tri.v[2] = old0.v[ino1_0];
    ref_tri.s2[0] = old0.s2[ino2_0];	ref_tri.s2[1] = itri2;			ref_tri.s2[2] = itri0;
    ////////////////
    if (old1.s2[ino2_0]>=0&&old1.s2[ino2_0]<(int)aTri.size()){
      assert(old0.r2[ino2_0] < 3);
      const unsigned int* rel = relTriTri[old0.r2[ino2_0]];
      ref_tri.r2[0] = noel2RelTriTri[rel[ino2_0]*3+rel[ino0_0]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old0.s2[ino2_0] < (int)aTri.size());
      aTri[old0.s2[ino2_0]].s2[rel[ino2_0]] = itri1;
      aTri[old0.s2[ino2_0]].r2[rel[ino2_0]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri2];
    ////////////////
    ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old1.v[ino2_1];	ref_tri.v[2] = old1.v[ino0_1];
    ref_tri.s2[0] = old1.s2[ino1_1];	ref_tri.s2[1] = itri3;			ref_tri.s2[2] = itri1;
    ////////////////
    if (old1.s2[ino1_1]>=0&&old1.s2[ino1_1]<(int)aTri.size()){
      assert(old1.r2[ino1_1] < 3);
      const unsigned int* rel = relTriTri[old1.r2[ino1_1]];
      ref_tri.r2[0] = noel2RelTriTri[rel[ino1_1]*3+rel[ino2_1]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old1.s2[ino1_1] < (int)aTri.size());
      aTri[old1.s2[ino1_1]].s2[rel[ino1_1]] = itri2;
      aTri[old1.s2[ino1_1]].r2[rel[ino1_1]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri3];
    ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old1.v[ino0_1];	ref_tri.v[2] = old1.v[ino1_1];
    ref_tri.s2[0] = old1.s2[ino2_1];	ref_tri.s2[1] = itri0;			ref_tri.s2[2] = itri2;
    if (old1.s2[ino2_1]>=0&&old1.s2[ino2_1]<(int)aTri.size()){
      assert(old1.r2[ino2_1] < 3);
      const unsigned int* rel = relTriTri[old1.r2[ino2_1]];
      ref_tri.r2[0] = noel2RelTriTri[rel[ino2_1]*3+rel[ino0_1]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old1.s2[ino2_1] < (int)aTri.size());
      aTri[old1.s2[ino2_1]].s2[rel[ino2_1]] = itri3;
      aTri[old1.s2[ino2_1]].r2[rel[ino2_1]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  return true;
}

template <typename TYPE>
bool InsertPoint_Elem
(const int ipo_ins,
const int itri_ins,
std::vector<CEPo<TYPE> >& aPo,
std::vector<ETri>& aTri)
{
  assert(itri_ins < (int)aTri.size());
  assert(ipo_ins < (int)aPo.size());
  
  // (node index opp to 0)*3+(node index opp to 1) -> relation index
  const int noel2RelTriTri[9] = {
    -1,	// 0 00
    -1,	// 1 01
    0,	// 2 02
    2, // 3 10
    -1, // 4 11
    -1,	// 5 12
    -1,	// 6 20
    1,	// 7 21
    -1, // 8 22
  };

  const int itri0 = itri_ins;
  const int itri1 = (int)aTri.size();
  const int itri2 = (int)aTri.size()+1;

  aTri.resize(aTri.size()+2);
  const ETri old_tri = aTri[itri_ins];

  aPo[ipo_ins].e = itri0;			  aPo[ipo_ins].d = 0;
  aPo[old_tri.v[0]].e = itri1;	aPo[old_tri.v[0]].d = 2;
  aPo[old_tri.v[1]].e = itri2;	aPo[old_tri.v[1]].d = 2;
  aPo[old_tri.v[2]].e = itri0;	aPo[old_tri.v[2]].d = 2;

  {
    ETri& ref_tri = aTri[itri0];
    ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old_tri.v[1];	ref_tri.v[2] = old_tri.v[2];
    ref_tri.s2[0] = old_tri.s2[0];	ref_tri.s2[1] = itri1;			ref_tri.s2[2] = itri2;
    if (old_tri.s2[0]>=0&&old_tri.s2[0]<(int)aTri.size()){
      assert(old_tri.r2[0] < 3);
      const unsigned int* rel = &relTriTri[old_tri.r2[0]][0];
      ref_tri.r2[0] = noel2RelTriTri[rel[0]*3+rel[1]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old_tri.s2[0] < (int)aTri.size());
      aTri[old_tri.s2[0]].s2[rel[0]] = itri0;
      aTri[old_tri.s2[0]].r2[rel[0]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri1];
    ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old_tri.v[2];	ref_tri.v[2] = old_tri.v[0];
    ref_tri.s2[0] = old_tri.s2[1];	ref_tri.s2[1] = itri2;			ref_tri.s2[2] = itri0;
    if (old_tri.s2[1]>=0&&old_tri.s2[1]<(int)aTri.size()){
      assert(old_tri.r2[1] < 3);
      const unsigned int* rel = &relTriTri[old_tri.r2[1]][0];
      ref_tri.r2[0] = noel2RelTriTri[rel[1]*3+rel[2]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old_tri.s2[1] < (int)aTri.size());
      aTri[old_tri.s2[1]].s2[rel[1]] = itri1;
      aTri[old_tri.s2[1]].r2[rel[1]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri2];
    ref_tri.v[0] = ipo_ins; 		ref_tri.v[1] = old_tri.v[0];	ref_tri.v[2] = old_tri.v[1];
    ref_tri.s2[0] = old_tri.s2[2];	ref_tri.s2[1] = itri0;			ref_tri.s2[2] = itri1;
    if (old_tri.s2[2]>=0&&old_tri.s2[2]<(int)aTri.size()){
      assert(old_tri.r2[2] < 3);
      const unsigned int* rel = &relTriTri[old_tri.r2[2]][0];
      ref_tri.r2[0] = noel2RelTriTri[rel[2]*3+rel[0]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old_tri.s2[2] < (int)aTri.size());
      aTri[old_tri.s2[2]].s2[rel[2]] = itri2;
      aTri[old_tri.s2[2]].r2[rel[2]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  return true;
}

template <typename TYPE>
int InsertPoint_Mesh
(const int itri0,
 double& r0,
 double& r1,
 std::vector<CEPo<TYPE> >& aPo3D,
 std::vector<ETri>& aSTri)
{
  if (itri0==-1) return -1;
  CVector3 pos,norm;
  {
    const int i0 = aSTri[itri0].v[0];
    const int i1 = aSTri[itri0].v[1];
    const int i2 = aSTri[itri0].v[2];
    const CVector3& p0 = aPo3D[i0].p;
    const CVector3& p1 = aPo3D[i1].p;
    const CVector3& p2 = aPo3D[i2].p;
    pos = r0*p0 + r1*p1 + (1-r0-r1)*p2;
    UnitNormal(norm, p0, p1, p2);
  }
  CEPo<void*> q(pos, norm);
  int ipo_ins = (int)aPo3D.size();
  aPo3D.push_back(q);
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

template <typename TYPE>
bool FlipEdge(int itri0, int ied0,
  std::vector<CEPo<TYPE> >& aPo, std::vector<ETri>& aTri)
{
  assert(itri0 < (int)aTri.size());
  assert(ied0 < 3);
  assert(aTri[itri0].s2[ied0]>=0&&aTri[itri0].s2[ied0]<(int)aTri.size());

  const int itri1 = aTri[itri0].s2[ied0];
  const int ied1 = relTriTri[aTri[itri0].r2[ied0]][ied0];
  assert(itri1 < (int)aTri.size());
  assert(ied1 < 3);
  assert(aTri[itri1].s2[ied1]>=0&&aTri[itri0].s2[ied0]<(int)aTri.size());
  
  // (node index opp to 0)*3+(node index opp to 1) -> relation index
  const int noel2RelTriTri[9] = {
    -1,	// 0 00
    -1,	// 1 01
    0,	// 2 02
    2, // 3 10
    -1, // 4 11
    -1,	// 5 12
    -1,	// 6 20
    1,	// 7 21
    -1, // 8 22
  };

  //	std::cout << itri0 << "-" << ied0 << "    " << itri1 << "-" << ied1 << std::endl;

  ETri old0 = aTri[itri0];
  ETri old1 = aTri[itri1];

  const int no0_0 = ied0;
  const int no1_0 = (ied0+1)%3; //noelTriEdge[ied0][0];
  const int no2_0 = (ied0+2)%3; //noelTriEdge[ied0][1];

  const int no0_1 = ied1;
  const int no1_1 = (ied1+1)%3; //noelTriEdge[ied1][0];
  const int no2_1 = (ied1+2)%3; //noelTriEdge[ied1][1];
  
  if(old0.s2[no1_0]==old1.s2[no2_1]) return false; // inverted mesh
  if(old0.s2[no2_0]==old1.s2[no1_1]) return false; // inverted mesh

  assert(old0.v[no1_0]==old1.v[no2_1]);
  assert(old0.v[no2_0]==old1.v[no1_1]);

  aPo[old0.v[no1_0]].e = itri0;	aPo[old0.v[no1_0]].d = 0;
  aPo[old0.v[no0_0]].e = itri0;	aPo[old0.v[no0_0]].d = 2;
  aPo[old1.v[no1_1]].e = itri1;	aPo[old1.v[no1_1]].d = 0;
  aPo[old1.v[no0_1]].e = itri1;	aPo[old1.v[no0_1]].d = 2;

  {
    ETri& ref_tri = aTri[itri0];
    ////////////////
    ref_tri.v[0] = old0.v[no1_0];	ref_tri.v[1] = old1.v[no0_1];	ref_tri.v[2] = old0.v[no0_0];
    ref_tri.s2[0] = itri1;			ref_tri.s2[1] = old0.s2[no2_0];	ref_tri.s2[2] = old1.s2[no1_1];
    ////////////////
    ref_tri.r2[0] = 0;
    if (old0.s2[no2_0]>=0){
      assert(old0.r2[no2_0] < 3);
      const unsigned int* rel = relTriTri[old0.r2[no2_0]];
      assert(old0.s2[no2_0] < (int)aTri.size());
      assert(old0.s2[no2_0]!=itri0);
      assert(old0.s2[no2_0]!=itri1);
      ref_tri.r2[1] = noel2RelTriTri[rel[no1_0]*3+rel[no2_0]];
      assert(ref_tri.r2[1]>=0&&ref_tri.r2[1] < 3);
      aTri[old0.s2[no2_0]].s2[rel[no2_0]] = itri0;
      aTri[old0.s2[no2_0]].r2[rel[no2_0]] = invRelTriTri[ref_tri.r2[1]];
    }
    if (old1.s2[no1_1]>=0){
      assert(old1.r2[no1_1] < 3);
      const unsigned int* rel = relTriTri[old1.r2[no1_1]];
      assert(old1.s2[no1_1] < (int)aTri.size());
      ref_tri.r2[2] = noel2RelTriTri[rel[no2_1]*3+rel[no0_1]];
      assert(ref_tri.r2[2]>=0&&ref_tri.r2[2] < 3);
      aTri[old1.s2[no1_1]].s2[rel[no1_1]] = itri0;
      aTri[old1.s2[no1_1]].r2[rel[no1_1]] = invRelTriTri[ref_tri.r2[2]];
    }
  }

  {
    ETri& ref_tri = aTri[itri1];
    ////////////////
    ref_tri.v[0] = old1.v[no1_1];	ref_tri.v[1] = old0.v[no0_0];	ref_tri.v[2] = old1.v[no0_1];
    ref_tri.s2[0] = itri0;			ref_tri.s2[1] = old1.s2[no2_1];	ref_tri.s2[2] = old0.s2[no1_0];
    ////////////////
    ref_tri.r2[0] = 0;
    if (old1.s2[no2_1]>=0){
      assert(old1.r2[no2_1] < 3);
      const unsigned int* rel = relTriTri[old1.r2[no2_1]];
      assert(old1.s2[no2_1] < (int)aTri.size());
      ref_tri.r2[1] = noel2RelTriTri[rel[no1_1]*3+rel[no2_1]];
      assert(ref_tri.r2[1]>=0&&ref_tri.r2[1] < 3);
      aTri[old1.s2[no2_1]].s2[rel[no2_1]] = itri1;
      aTri[old1.s2[no2_1]].r2[rel[no2_1]] = invRelTriTri[ref_tri.r2[1]];
    }
    if (old0.s2[no1_0]>=0){
      assert(old0.r2[no1_0] < 3);
      const unsigned int* rel = relTriTri[old0.r2[no1_0]];
      assert(old0.s2[no1_0] < (int)aTri.size());
      ref_tri.r2[2] = noel2RelTriTri[rel[no2_0]*3+rel[no0_0]];
      assert(ref_tri.r2[2]>=0&&ref_tri.r2[2] < 3);
      aTri[old0.s2[no1_0]].s2[rel[no1_0]] = itri1;
      aTri[old0.s2[no1_0]].r2[rel[no1_0]] = invRelTriTri[ref_tri.r2[2]];
    }
  }

  return true;
}

template <typename TYPE>
void MoveCCW(int& itri_cur, int& inotri_cur, bool& flag_is_wall,
//             const int itri0, const int ipo0,
             std::vector<CEPo<TYPE> >& aPo, std::vector<ETri>& aTri)
{
  const int inotri1 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
  if (aTri[itri_cur].s2[inotri1]==-1){ flag_is_wall = true; return; }
  flag_is_wall = false;
  const int itri_nex = aTri[itri_cur].s2[inotri1];
  const unsigned int* rel_nex = relTriTri[aTri[itri_cur].r2[inotri1]];
  const int inotri_nex = rel_nex[inotri_cur];
  itri_cur = itri_nex;
  inotri_cur = inotri_nex;
}

template <typename TYPE>
bool DelaunayAroundPoint(int ipo0,
  std::vector<CEPo<TYPE> >& aPo, std::vector<ETri>& aTri)
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
      if (DetDelaunay(
        aPo[aTri[itri_cur].v[0]].p,
        aPo[aTri[itri_cur].v[1]].p,
        aPo[aTri[itri_cur].v[2]].p,
        aPo[ipo_dia].p)==0)
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
    {	// next element
      const int inotri1 = indexRot3[1][inotri_cur];
      if (aTri[itri_cur].s2[inotri1]==-1){
        flag_is_wall = true;
        break;
      }
      const int itri_nex = aTri[itri_cur].s2[inotri1];
      const unsigned int* rel_nex = relTriTri[aTri[itri_cur].r2[inotri1]];
      const int inotri_nex = rel_nex[inotri_cur];
      assert(aTri[itri_nex].v[inotri_nex]==ipo0);
      if (itri_nex==itri0) break;	// finish if we reach starting elemnt
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
      if (DetDelaunay(
        aPo[aTri[itri_cur].v[0]].p,
        aPo[aTri[itri_cur].v[1]].p,
        aPo[aTri[itri_cur].v[2]].p,
        aPo[ipo_dia].p)==0)	// Delaunay condition is not satisfiled
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
      assert(itri_nex!=itri0);	// finsih if reach starting elemnet
      itri_cur = itri_nex;
      inotri_cur = inotri_nex;
    }
  }
  return true;
}

template <typename TYPE>
bool DeleteTri
(unsigned int itri_to,
std::vector<CEPo<TYPE> >& aPo,
std::vector<ETri>& aTri)
{
  if (itri_to>=aTri.size()) return true;
  {
    assert(aTri[itri_to].s2[0]==-1);
    assert(aTri[itri_to].s2[1]==-1);
    assert(aTri[itri_to].s2[2]==-1);
  }
  const unsigned int itri_from = (int)aTri.size()-1;
  if (itri_to==itri_from){
    aTri.resize(aTri.size()-1);
    return true;
  }
  aTri[itri_to] = aTri[itri_from];
  aTri.resize(aTri.size()-1);
  for (int iedtri = 0; iedtri<3; iedtri++){
    if (aTri[itri_to].s2[iedtri]==-1) continue;
    const unsigned int itri_adj = aTri[itri_to].s2[iedtri];
    assert(itri_adj < aTri.size());
    const unsigned int* rel = relTriTri[(int)aTri[itri_to].r2[iedtri]];
    const unsigned int iedtri_adj = rel[iedtri];
    assert(aTri[itri_adj].s2[iedtri_adj]==itri_from);
    aTri[itri_adj].s2[iedtri_adj] = itri_to;
  }
  for (unsigned int inotri = 0; inotri<3; inotri++){
    const unsigned int ipo0 = aTri[itri_to].v[inotri];
    aPo[ipo0].e = itri_to;
    aPo[ipo0].d = inotri;
  }
  return true;
}

template <typename TYPE>
bool Collapse_ElemEdge
(const int itri_del, const int ied_del,
std::vector<CEPo<TYPE> >& aPo, std::vector<ETri>& aTri)
{
  assert(itri_del < (int)aTri.size());
  if (aTri[itri_del].s2[ied_del]==-1){
    std::cout<<"Error!-->Not Implemented: Mesh with hole"<<std::endl;
    assert(0);
  }

  const int itri_adj = aTri[itri_del].s2[ied_del];
  const int ied_adj = (int)relTriTri[(int)aTri[itri_del].r2[ied_del]][ied_del];
  assert(itri_adj < (int)aTri.size());
  assert(ied_adj < 3);
  assert(aTri[itri_adj].s2[ied_adj]==itri_del);

  const int itri0 = itri_del;
  const int itri1 = itri_adj;
  const int itri2 = aTri[itri0].s2[(ied_del+1)%3];
  const int itri3 = aTri[itri0].s2[(ied_del+2)%3];
  const int itri4 = aTri[itri1].s2[(ied_adj+1)%3];
  const int itri5 = aTri[itri1].s2[(ied_adj+2)%3];

  const int ino0_0 = ied_del;
  const int ino1_0 = (ino0_0+1)%3;
  const int ino2_0 = (ino0_0+2)%3;

  const int ino0_1 = ied_adj;
  const int ino1_1 = (ino0_1+1)%3;
  const int ino2_1 = (ino0_1+2)%3;

  const unsigned int ino0_2 = relTriTri[(int)aTri[itri0].r2[ino1_0]][ino1_0];
  const int ino1_2 = (ino0_2+1)%3;
  const int ino2_2 = (ino0_2+2)%3;

  const unsigned ino0_3 = relTriTri[(int)aTri[itri0].r2[ino2_0]][ino2_0];
  const int ino1_3 = (ino0_3+1)%3;

  const unsigned ino0_4 = relTriTri[(int)aTri[itri1].r2[ino1_1]][ino1_1];
  const int ino1_4 = (ino0_4+1)%3;
  const int ino2_4 = (ino0_4+2)%3;

  const unsigned ino0_5 = relTriTri[(int)aTri[itri1].r2[ino2_1]][ino2_1];
  const int ino1_5 = (ino0_5+1)%3;

  if (aTri[itri2].s2[ino2_2]==itri3) return false;
  if (aTri[itri3].s2[ino1_3]==itri2) return false;
  if (aTri[itri4].s2[ino2_4]==itri5) return false;
  if (aTri[itri5].s2[ino1_5]==itri4) return false;
  if (itri2==itri5 && itri3==itri4) return false;

  const ETri old0 = aTri[itri0];
  const ETri old1 = aTri[itri1];

  int ipo0 = old0.v[ino0_0];
  int ipo1 = old0.v[ino1_0];
  int ipo2 = old1.v[ino0_1];
  int ipo3 = old1.v[ino1_1];  // to be delete

  {
    std::vector<int> ring1;
    { // set triangle index from point 0 to point 1
      int jtri = itri5;
      int jnoel_c = ino1_5;
      int jnoel_b = (jnoel_c+1)%3; // noelTriEdge[jnoel_c][0];
      for (;;){
        assert(jtri < (int)aTri.size());
        assert(jnoel_c < 3);
        assert(aTri[jtri].v[jnoel_c]==ipo3);
        {
          int jpo = aTri[jtri].v[(jnoel_c+2)%3];
          ring1.push_back(jpo);
        }
        assert(aTri[jtri].s2[jnoel_b]>=0&&aTri[jtri].s2[jnoel_b]<(int)aTri.size());
        int ktri = aTri[jtri].s2[jnoel_b];
        const int rel01 = aTri[jtri].r2[jnoel_b];
        const int knoel_c = relTriTri[rel01][jnoel_c];
        const int knoel_b = relTriTri[rel01][(jnoel_c+2)%3];
        assert(itri1 < (int)aTri.size());
        assert(aTri[ktri].s2[relTriTri[rel01][jnoel_b]]==jtri);
        if (ktri==itri2) break;
        jtri = ktri;
        jnoel_c = knoel_c;
        jnoel_b = knoel_b;
      }
    }
    std::vector<int> ring2;
    { // set triangle index from point 0 to point 1
      int jtri = itri3;
      int jnoel_c = ino1_3;
      int jnoel_b = (jnoel_c+1)%3; // noelTriEdge[jnoel_c][0];
      for (;;){
        assert(jtri < (int)aTri.size());
        assert(jnoel_c < 3);
        assert(aTri[jtri].v[jnoel_c]==ipo1);
        {
          int jpo = aTri[jtri].v[(jnoel_c+2)%3];
          ring2.push_back(jpo);
        }
        assert(aTri[jtri].s2[jnoel_b]>=0&&aTri[jtri].s2[jnoel_b]<(int)aTri.size());
        int ktri = aTri[jtri].s2[jnoel_b];
        const int rel01 = aTri[jtri].r2[jnoel_b];
        const int knoel_c = relTriTri[rel01][jnoel_c];
        const int knoel_b = relTriTri[rel01][(jnoel_c+2)%3];
        assert(itri1 < (int)aTri.size());
        assert(aTri[ktri].s2[relTriTri[rel01][jnoel_b]]==jtri);
        if (ktri==itri4) break;
        jtri = ktri;
        jnoel_c = knoel_c;
        jnoel_b = knoel_b;
      }
    }
    sort(ring1.begin(), ring1.end());
    sort(ring2.begin(), ring2.end());
    std::vector<int> insc(ring1.size());
    std::vector<int>::iterator it = set_intersection(ring1.begin(), ring1.end(), ring2.begin(), ring2.end(), insc.begin());
    if (it!=insc.begin()){ return  false; }
  }


  std::cout<<std::endl;
  std::cout<<"stt"<<std::endl;
  std::cout<<"tris : "<<itri0<<" "<<itri1<<" "<<itri2<<" "<<itri3<<" "<<itri4<<" "<<itri5<<std::endl;
  std::cout<<"vtxs : "<<ipo0<<" "<<ipo1<<" "<<ipo2<<" "<<ipo3<<std::endl;

  assert(old0.v[ino1_0]==old1.v[ino2_1]);
  assert(old0.v[ino2_0]==old1.v[ino1_1]);
  assert(old0.s2[ino0_0]==itri1);
  assert(old1.s2[ino0_1]==itri0);

  {
    aPo[ipo0].e = itri2;	aPo[ipo0].d = ino1_2;
    aPo[ipo2].e = itri4;	aPo[ipo2].d = ino1_4;
    aPo[ipo1].e = itri3;	aPo[ipo1].d = ino1_3;
    aPo[ipo3].e = -1;
  }
  
  // (edge index)*3+(opp edge index) -> relation index
  const int ed2RelTriTri[9] = {
    0,	// 0 00
    2,	// 1 01
    1,	// 2 02
    2,	// 3 10
    1,	// 4 11
    0,	// 5 12
    1,	// 6 20
    0,	// 7 21
    2,	// 8 22
  };

  { // change itri2
    ETri& tri = aTri[itri2];
    //    tri.g2[ino0_2] = old0.g2[ino2_0];
    tri.s2[ino0_2] = old0.s2[ino2_0];
    if (old0.s2[ino2_0]>=0&&old0.s2[ino2_0]<(int)aTri.size()){
      assert(old0.r2[ino2_0] < 3);
      assert(old0.s2[ino2_0] < (int)aTri.size());
      tri.r2[ino0_2] = ed2RelTriTri[ino0_2*3+ino0_3];
      assert(tri.r2[ino0_2]>=0&&tri.r2[ino0_2] < 3);
      aTri[itri3].s2[ino0_3] = itri2;
      aTri[itri3].r2[ino0_3] = invRelTriTri[tri.r2[ino0_2]];
    }
  }

  { // change itri3
    ETri& tri = aTri[itri3];
    //    tri.g2[ino0_3] = old0.g2[ino1_0];
    tri.s2[ino0_3] = old0.s2[ino1_0];
    if (old0.s2[ino1_0]>=0&&old0.s2[ino1_0]<(int)aTri.size()){
      //		if( old0.g2[ino1_0] == -2 || old0.g2[ino1_0] == -3 ){
      assert(old0.r2[ino1_0] < 3);
      assert(old0.s2[ino1_0] < (int)aTri.size());
      tri.r2[ino0_3] = ed2RelTriTri[ino0_3*3+ino0_2];
      assert(tri.r2[ino0_3]>=0&&tri.r2[ino0_3] < 3);
      aTri[itri2].s2[ino0_2] = itri3;
      aTri[itri2].r2[ino0_2] = invRelTriTri[tri.r2[ino0_3]];
    }
  }

  { // change itri4
    ETri& tri = aTri[itri4];
    //    tri.g2[ino0_4] = old1.g2[ino2_1];
    tri.s2[ino0_4] = old1.s2[ino2_1];
    //		if( old1.g2[ino2_1] == -2 || old1.g2[ino2_1] == -3 ){
    if (old1.s2[ino2_1]>=0&&old1.s2[ino2_1] < (int)aTri.size()){
      assert(old1.r2[ino2_1] < 3);
      assert(old1.s2[ino2_1] < (int)aTri.size());
      tri.r2[ino0_4] = ed2RelTriTri[ino0_4*3+ino0_5];
      assert(tri.r2[ino0_4]>=0&&tri.r2[ino0_4] < 3);
      aTri[itri5].s2[ino0_5] = itri4;
      aTri[itri5].r2[ino0_5] = invRelTriTri[tri.r2[ino0_4]];
      assert(relTriTri[(int)aTri[itri4].r2[ino0_4]][ino0_4]==ino0_5);
      assert(relTriTri[(int)aTri[itri5].r2[ino0_5]][ino0_5]==ino0_4);
    }
  }

  { // change itri5
    ETri& tri = aTri[itri5];
    //    tri.g2[ino0_5] = old1.g2[ino1_1];
    tri.s2[ino0_5] = old1.s2[ino1_1];
    //		if( old1.g2[ino1_1] == -2 || old1.g2[ino1_1] == -3 ){
    if (old1.s2[ino1_1]>=0&&old1.s2[ino1_1] < (int)aTri.size()){
      assert(old1.r2[ino1_1] < 3);
      assert(old1.s2[ino1_1] < (int)aTri.size());
      tri.r2[ino0_5] = ed2RelTriTri[ino0_5*3+ino0_4];
      assert(tri.r2[ino0_5]>=0&&tri.r2[ino0_5] < 3);
      aTri[itri4].s2[ino0_4] = itri5;
      aTri[itri4].r2[ino0_4] = invRelTriTri[tri.r2[ino0_5]];
      assert(relTriTri[aTri[itri5].r2[ino0_5]][ino0_5]==ino0_4);
      assert(relTriTri[aTri[itri4].r2[ino0_4]][ino0_4]==ino0_5);
    }
  }

  { // set triangle index from point 0 to point 1
    int jtri = itri5;
    int jnoel_c = ino1_5;
    int jnoel_b = (jnoel_c+1)%3; // noelTriEdge[jnoel_c][0];
    for (;;){
      assert(jtri < (int)aTri.size());
      assert(jnoel_c < 3);
      assert(aTri[jtri].v[jnoel_c]==ipo3);
      //      std::cout << " fan : " << jtri << "     " << aTri[jtri].v[0] << " " << aTri[jtri].v[1] << " " << aTri[jtri].v[2] << std::endl;
      aTri[jtri].v[jnoel_c] = ipo1;
      //      assert( aTri[jtri].g2[jnoel_b] == -2 );
      assert(aTri[jtri].s2[jnoel_b]>=0&&aTri[jtri].s2[jnoel_b]<(int)aTri.size());
      int ktri = aTri[jtri].s2[jnoel_b];
      const int rel01 = aTri[jtri].r2[jnoel_b];
      const int knoel_c = relTriTri[rel01][jnoel_c];
      const int knoel_b = relTriTri[rel01][(jnoel_c+2)%3];
      assert(itri1 < (int)aTri.size());
      assert(aTri[ktri].s2[relTriTri[rel01][jnoel_b]]==jtri);
      if (ktri==itri3||ktri==itri4) break;
      jtri = ktri;
      jnoel_c = knoel_c;
      jnoel_b = knoel_b;
    }
  }

  {	// isolate two triangles to be deleted
    aTri[itri0].s2[0] = -1;  aTri[itri0].s2[1] = -1;  aTri[itri0].s2[2] = -1;
    aTri[itri1].s2[0] = -1;  aTri[itri1].s2[1] = -1;  aTri[itri1].s2[2] = -1;
    const int itri_1st = (itri0 > itri1) ? itri0 : itri1;
    const int itri_2nd = (itri0 < itri1) ? itri0 : itri1;
    DeleteTri(itri_1st, aPo, aTri);
    DeleteTri(itri_2nd, aPo, aTri);
  }
  return true;
}

// information look-up
template <typename TYPE>
void GetTriAryAroundPoint
(int ipoin,
 std::vector< std::pair<unsigned int,unsigned int> >& aTriSurPo,
 const std::vector<CEPo<TYPE> >& aPo, const std::vector<ETri>& aTri)
{
  const unsigned int itri_ini = aPo[ipoin].e;
  const unsigned int inoel_c_ini = aPo[ipoin].d;
  assert(itri_ini < aTri.size());
  assert(inoel_c_ini < 3);
  assert(aTri[itri_ini].v[inoel_c_ini]==ipoin);
  unsigned int itri0 = itri_ini;
  unsigned int inoel_c0 = inoel_c_ini;
  unsigned int inoel_b0 = (inoel_c0+1)%3;
  for (;;){
    assert(itri0 < aTri.size());
    assert(inoel_c0 < 3);
    assert(aTri[itri0].v[inoel_c0]==ipoin);
    aTriSurPo.push_back(std::make_pair(itri0, inoel_c0));
    /////
    //    std::cout << ipoin << " " << itri0 << " " << inoel_b0 << " " << aTri.size() << " " << aTri[itri0].s2[inoel_b0] << std::endl;
    if (aTri[itri0].s2[inoel_b0]==-1){ break; }
    assert(aTri[itri0].s2[inoel_b0]>=0&&aTri[itri0].s2[inoel_b0] < (int)aTri.size());
    unsigned int itri1 = aTri[itri0].s2[inoel_b0];
    const unsigned int rel01 = aTri[itri0].r2[inoel_b0];
    const unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
    const unsigned int inoel_b1 = relTriTri[rel01][(inoel_c0+2)%3];
    assert(itri1 < aTri.size());
    assert(aTri[itri1].s2[relTriTri[rel01][inoel_b0]]==itri0);
    assert(aTri[itri1].v[inoel_c1]==ipoin);
    if (itri1==itri_ini) return;
    itri0 = itri1;
    inoel_c0 = inoel_c1;
    inoel_b0 = inoel_b1;
  }
}

template <typename TYPE>
bool FindRayTriangleMeshIntersections
(const CVector3 &line0,
 const CVector3 &line1,
 const std::vector<ETri>& aTri,
 const std::vector<CEPo<TYPE> > &aPoint3D,
 std::vector<CVector3> &intersectionPoints);

template <typename TYPE>
bool FindRayTriangleMeshIntersectionClosestToPoint
(const CVector3 &line0,
 const CVector3 &line1,
 const std::vector<ETri>& aTri,
 const std::vector<CEPo<TYPE> > &aPoint3D,
 const CVector3 &targetPoint,
 CVector3 &intersectionPoint);


/////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TYPE>
void InitializeMesh
(std::vector<CEPo<TYPE> >& aPo3D,
std::vector<ETri>& aSTri,
////
const std::vector<double>& aXYZ,
const std::vector<int>& aTri)
{
  aPo3D.clear();
  aSTri.clear();
  /////
  aPo3D.resize(aXYZ.size()/3);
//  std::cout << aPo3D.size() << std::endl;
  for (int ipo = 0; ipo<(int)aPo3D.size(); ++ipo){
    aPo3D[ipo].p.x = aXYZ[ipo*3+0];
    aPo3D[ipo].p.y = aXYZ[ipo*3+1];
    aPo3D[ipo].p.z = aXYZ[ipo*3+2];
    aPo3D[ipo].e = -1; // for unreffered point
    aPo3D[ipo].d = 0;
  }
  aSTri.resize(aTri.size()/3);
  for (int itri = 0; itri<(int)aSTri.size(); itri++){
    aSTri[itri].v[0] = aTri[itri*3+0];
    aSTri[itri].v[1] = aTri[itri*3+1];
    aSTri[itri].v[2] = aTri[itri*3+2];
  }
  for (int itri = 0; itri<(int)aSTri.size(); itri++){
    unsigned int i1 = aSTri[itri].v[0];
    unsigned int i2 = aSTri[itri].v[1];
    unsigned int i3 = aSTri[itri].v[2];
    aPo3D[i1].e = itri; aPo3D[i1].d = 0;
    aPo3D[i2].e = itri; aPo3D[i2].d = 1;
    aPo3D[i3].e = itri; aPo3D[i3].d = 2;
  }
  {
    unsigned int* elsup_ind = new unsigned int[aPo3D.size()+1];
    unsigned int nelsup;
    unsigned int* elsup;
    MakePointSurTri(aSTri, (int)aPo3D.size(), elsup_ind, nelsup, elsup);
    MakeInnerRelationTri(aSTri, (int)aPo3D.size(), elsup_ind, nelsup, elsup);
    delete[] elsup_ind;
    delete[] elsup;
  }
  makeNormal(aPo3D, aSTri);
  //  CheckTri(aSTri);
  //  CheckTri(aPo3D, aSTri);
}


template <typename TYPE>
int pickTriangle
(CVector3& p,
const CVector3& org, const CVector3& dir,
int itri_start, // starting triangle
const std::vector<CEPo<TYPE> >& aPo3D,
const std::vector<ETri>& aSTri)
{
  int itri1 = itri_start;
  for (int itr = 0; itr<50; itr++){
    if (itri1==-1) return -1;
    int ip0 = aSTri[itri1].v[0];
    int ip1 = aSTri[itri1].v[1];
    int ip2 = aSTri[itri1].v[2];
    const CVector3& tp0 = aPo3D[ip0].p;
    const CVector3& p1 = aPo3D[ip1].p;
    const CVector3& p2 = aPo3D[ip2].p;
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


template <typename TYPE>
bool pickMesh(CVector3& p,
              int& itriOut,
              double& r0Out,
              double& r1Out,
              ////
              const CVector3& org,
              const CVector3& dir,
              int itri_start, // starting triangle
              const std::vector<CEPo<TYPE> >& aPo3D,
              const std::vector<ETri>& aSTri)
{
  int itri1 = itri_start;
  for (int itr = 0; itr<50; itr++){
    if (itri1==-1) return false;
    int ip0 = aSTri[itri1].v[0];
    int ip1 = aSTri[itri1].v[1];
    int ip2 = aSTri[itri1].v[2];
    const CVector3& tp0 = aPo3D[ip0].p;
    const CVector3& p1 = aPo3D[ip1].p;
    const CVector3& p2 = aPo3D[ip2].p;
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

void extractHoles
(std::vector< std::vector<int> >& aIndP_Hole,
 const int npo,
 const std::vector<ETri>& aETri);

//////////////////////////////////////////////////////////////////////////////////////////////////

class CMeshDensity
{
public:
  virtual double edgeLengthRatio(double px, double py) const = 0;
};

int delaunay_triangulation2(std::vector<int>& aTri_out,		// out
                            std::vector<double>& aXY_out, // out
                            std::vector<int>& aPtrVtxInd, // out
                            std::vector<int>& aVtxInd, // out
                            ////
                            bool is_add_point_boundary,
                            const std::vector<int>& aIndXYs, // in
                            const std::vector<double>& aXY_in, // ind
                            const double max_edge_length, // ind
                            const CMeshDensity& mesh_density); // ind

// TODO: there should be three optoions for adding point on edge, 0:none, 1:only between input points, 2: resample everything
bool GenerateTesselation2(std::vector<int>& aTri_out, // out
                          std::vector<double>& aXY_out, // out
                          std::vector<int>& aPtrVtxInd,
                          std::vector<int>& aVtxInd,
                          ////
                          double elen,
                          const CMeshDensity& mesh_density,
                          bool is_uniform_resample_loop, // good for polyline curve in
                          const std::vector< std::vector<double> >& aVecAry0); // in



bool GenerateTesselation2(std::vector<int>& aTri_out, // out
                          std::vector<double>& aXY_out, // out
                          std::vector<int>& aPtrVtxInd,
                          std::vector<int>& aVtxInd,
                          ////
                          double elen,
                          bool is_uniform_resample_loop, // good for polyline curve
                          const std::vector< std::vector<double> >& aVecAry0); // in


#ifdef USE_OPENGL

#if defined(__APPLE__) && defined(__MACH__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#else
#include <windows.h>
#include <GL/gl.h>
#endif

template <typename TYPE>
void drawEditableMesh
(const std::vector< CEPo<TYPE> >& aPo3D,
 const std::vector<ETri>& aSTri )
{
  //  ::glPushAttrib(GL_ENABLE_BIT);
  
  ::glEnable(GL_LIGHTING);
  //  ::myGlColorDiffuse(CColor::Orange());
  ::glBegin(GL_TRIANGLES);
  for (int itri=0; itri<aSTri.size(); ++itri){
    const int i0 = aSTri[itri].v[0];
    const int i1 = aSTri[itri].v[1];
    const int i2 = aSTri[itri].v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
      continue;
    }
    {
      CVector3 n; UnitNormal(n, aPo3D[i0].p, aPo3D[i1].p, aPo3D[i2].p);
      ::glNormal3d(n.x,n.y,n.z);
    }
    {
      CVector3 p0 = aPo3D[i0].p;
      ::glVertex3d(p0.x,p0.y,p0.z);
    }
    {
      CVector3 p1 = aPo3D[i1].p;
      ::glVertex3d(p1.x,p1.y,p1.z);
    }
    {
      CVector3 p2 = aPo3D[i2].p;
      ::glVertex3d(p2.x,p2.y,p2.z);
    }
  }
  ::glEnd();
  
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (int itri = 0; itri<aSTri.size(); ++itri){
    const int i0 = aSTri[itri].v[0];
    const int i1 = aSTri[itri].v[1];
    const int i2 = aSTri[itri].v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
    }
    myGlVertex(aPo3D[i0].p);     myGlVertex(aPo3D[i1].p);
    myGlVertex(aPo3D[i1].p);     myGlVertex(aPo3D[i2].p);
    myGlVertex(aPo3D[i2].p);     myGlVertex(aPo3D[i0].p);
  }
  ::glEnd();
}

#endif


//////////////////////////////////////////////////////////////////////////////////////////////////


#endif // #endif SURFACE_MESH_H
