/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_DTET_V3_H
#define DFM2_DTET_V3_H

#include <math.h>
#include <vector>
#include <cassert>
#include <map>
#include <climits> // for UINT_MAX

#include "delfem2/dfm2_inline.h"
#include "delfem2/geosolidelm_v3.h"
#include "delfem2/geodelaunay3_v3.h"
#include "delfem2/vec3.h"

namespace delfem2 {

const double ILL_CRT = 1.0e+20;

/*! 
四面体のある面の番号
外からみて半時計周りになるように番号が並べられている．
*/
const int noelTetFace[4][3] = {
	{ 1, 2, 3 },
	{ 0, 3, 2 },
	{ 0, 1, 3 },
	{ 0, 2, 1 },
};

//! 四面体の隣接関係
static const unsigned int tetRel[12][4] = {
	{ 0, 1, 3, 2 }, //  0
	{ 0, 3, 2, 1 }, //  1 
	{ 0, 2, 1, 3 }, //  2

	{ 1, 0, 2, 3 }, //  3
	{ 1, 2, 3, 0 }, //  4
	{ 1, 3, 0, 2 }, //  5

	{ 2, 0, 3, 1 }, //  6
	{ 2, 3, 1, 0 }, //  7
	{ 2, 1, 0, 3 }, //  8

	{ 3, 0, 1, 2 }, //  9
	{ 3, 1, 2, 0 }, // 10
	{ 3, 2, 0, 1 }, // 11
};

//! 四面体の隣接関係の逆
const unsigned int invTetRel[12] = {
	 0, //  0
	 1, //  1
	 2, //  2
	 3, //  3
	 9, //  4
	 6, //  5
	 5, //  6 
	11, //  7
	 8, //  8
	 4, //  9
	10, // 10
	 7, // 11
};

//! 頂点０の対応関係aと頂点１の対応関係bにおいてa*4+bを引数として隣接関係を得る配列
const int noel2Rel[16] = {
	 -1, //  0 00
	  0, //  1 01
	  2, //  2 02
	  1, //  3 03
	  3, //  4 10
	 -1, //  5 11
	  4, //  6 12
	  5, //  7 13
	  6, //  8 20
	  8, //  9 21
	 -1, // 10 22
	  7, // 11 23
	  9, // 12 30
	 10, // 13 31
	 11, // 14 32
	 -1, // 15 33
};

//! 四面体中の無向辺の数
const unsigned int nSEdgeTet = 6;
//! 四面体中の無向辺の頂点(辞書式に並んでいる)
const unsigned int sEdge2Noel[6][2] = {
	{ 0, 1 },
	{ 0, 2 },
	{ 0, 3 },
	{ 1, 2 },
	{ 1, 3 },
	{ 2, 3 }
};
//! 四面体の無向辺の両側の面番号
const unsigned int sEdge2FaTet[6][2] = {
	{ 3, 2 },
	{ 1, 3 },
	{ 2, 1 },
	{ 3, 0 },
	{ 0, 2 },
	{ 1, 0 }
};

//! 四面体中の有向辺の数
const unsigned int nDEdgeTet = 12;
//! 四面体中の有向辺の頂点
const unsigned int dEdge2Noel[12][2] = {
	{ 0, 1 },	//  0 
	{ 0, 2 },	//  1
	{ 0, 3 },	//  2
	{ 1, 0 },	//  3
	{ 1, 2 },	//  4
	{ 1, 3 },	//  5
	{ 2, 0 },	//  6
	{ 2, 1 },	//  7
	{ 2, 3 },	//  8
	{ 3, 0 },	//  9
	{ 3, 1 },	// 10
	{ 3, 2 },	// 11
};
//! 四面体中の２頂点abを結ぶ有向辺の番号をa*4+bから得る配列
const int noel2DEdge[16] = {
	 -1, //  0 00
	  0, //  1 01	
	  1, //  2 02
	  2, //  3 03
	  3, //  4 10
	 -1, //  5 11
	  4, //  6 12
	  5, //  7 13
	  6, //  8 20
	  7, //  9 21
	 -1, // 10 22
	  8, // 11 23
	  9, // 12 30
	 10, // 13 31
	 11, // 14 32
	 -1, // 15 33
};
//! 無向辺から有向辺へのマッピング配列
const unsigned int sEdge2DEdge[6] = {
	0,
	1,
	2,
	4,
	5,
	8
};

////////////////////////////////////////////////////////////////
/*!
四面体のある面の番号
外からみて半時計周りになるように番号が並べられている．
*/
const unsigned int noelHexFace[6][4] = {
	{ 0, 3, 2, 1 },
	{ 0, 1, 5, 4 },
	{ 1, 2, 6, 5 }, 
	{ 2, 3, 7, 6 },
	{ 3, 0, 4, 7 },
	{ 4, 5, 6, 7 }
};

// -------------------------------------------

/**
 * @brief element around edge
 */
class ElemAroundEdge{
public:
	ElemAroundEdge(){
		e.reserve(32);
		n.reserve(32);
	}
	void clear(){
		e.clear();
		n.clear();
	}
  int size() const {
		return (int)e.size();
	}
public:
  std::vector< std::pair<unsigned int,unsigned int> > e;
  std::vector< unsigned int > n;
  unsigned int nod;
  unsigned int nou;
  bool is_inner;
  
};

/**
 * @brief point around element
 */
class ElemAroundPoint{
public:
	std::map<int,int> e;
	bool is_inner;
public:
	void clear(){
		e.clear();
	}
	int size() const {
		return (int)e.size();
	}
};

// -------------------------

/*
class CFlipCrtPrePosPtn
{
public:
	double pre;
	double pos;
	int ptn;
	bool ok_flg;
public:
	CFlipCrtPrePosPtn(const CFlipCrtPrePosPtn& rhs){
		this->pre = rhs.pre;
		this->pos = rhs.pos;
		this->ptn = rhs.ptn;
		this->ok_flg = rhs.ok_flg;
	}
	CFlipCrtPrePosPtn( double pre, double pos, int ptn ){
		this->pre = pre;	this->pos = pos;	this->ptn = ptn;
		if( pos > pre || ptn < 0 || pos > ILL_CRT*0.9 ) ok_flg = false;
		else ok_flg = true;
	}
	CFlipCrtPrePosPtn(){
		this->pre = -1.0;	this->pos = ILL_CRT;	this->ptn = -1;
		ok_flg = false;
	}
	void SetData( double pre, double pos, int ptn ){
		this->pre = pre;	this->pos = pos;	this->ptn = ptn;
		if( pos > pre || ptn < 0 || pos > ILL_CRT*0.9 ) ok_flg = false;
		else ok_flg = true;
	}
	bool operator < (const CFlipCrtPrePosPtn& rhs){
		if(  this->ok_flg && !rhs.ok_flg ) return true;
		if( !this->ok_flg &&  rhs.ok_flg ) return false;
		if( fabs( this->pre - rhs.pre ) > 1.0e-8 ){
			return this->pre > rhs.pre;
		}
		if( fabs( this->pos - rhs.pos ) > 1.0e-8 ){
			return this->pos < rhs.pos;
		}
		return true;
	}
	bool operator > (const CFlipCrtPrePosPtn& rhs){
		return !( *this < rhs );
	}
};

bool operator < (const CFlipCrtPrePosPtn& lhs, const CFlipCrtPrePosPtn& rhs);
*/


/**
 * @brief 3D point class for dynamic tetrahedra mesh
 */
class CDynPointTet{
public:
	CDynPointTet(){}
	CDynPointTet( const CDynPointTet& rhs ) 
        : e(rhs.e), poel(rhs.poel), p(rhs.p)//,old_p(rhs.old_p){}
  {}
	CDynPointTet(double x, double y, double z) 
        : e(UINT_MAX), poel(0), p(x,y,z)//, old_p(-1){}
  {}
	bool operator < (const CDynPointTet& rhs){
		if( fabs( this->p.x - rhs.p.x ) > 1.0e-5 ){
			return this->p.x < rhs.p.x;
		}
		if( fabs( this->p.y - rhs.p.y ) > 1.0e-5 ){
			return this->p.y < rhs.p.y;
		}
		if( fabs( this->p.z - rhs.p.z ) > 1.0e-5 ){
			return this->p.z < rhs.p.z;
		}
		return false;
	}
public:
	unsigned int e;		//< element nubmer
	unsigned int poel;	// element point number
  //
  delfem2::CVec3d p;
//  int old_p;	//!< 変更前の節点番号
};

/**
 * @brief class of dynamic tetrahedra
 */
class CDynTet
{
public:
  inline bool isActive() const { return v[0] != UINT_MAX; }
public:
  /**
   * @brief index of vertex
   */
	unsigned int v[4];
  /**
   * @brief ndex of adjacent element
   */
	unsigned int s[4];
};


int GetRelationshipTet(const unsigned int* t0,
                       const unsigned int* t1);

bool IsInsideCircumSphere(
    const delfem2::CVec3d& p,
    const delfem2::CDynTet t,
    const delfem2::CVec3d& c,
    const std::vector<CDynPointTet>& aPo3D);

/*
//! ６面体要素構造体
struct SHex{
	unsigned int v[8];	//<! 頂点のIndex
	unsigned int s[6];	//<! 隣接要素のIndex
	int g[6];
	unsigned int f[6];
	int old_h;
};

//! ３次元三角形要素構造体
struct STri3D{
	unsigned int v[3];	//<! 頂点のIndex
	////////////////
  unsigned int s2[3];	//<! 隣接要素のIndex
	int g2[3];
	unsigned int r2[3];
	////////////////
	int sf[2];	//<! 隣接要素のIndex
	int gf[2];
	unsigned int df[2];
};

//! ３次元４角形要素構造体
struct SQuad3D{
	unsigned int v[4];	//<! 頂点のIndex
	////////////////
	int se[4];	//<! 隣接要素のIndex
	int ge[4];
	unsigned int de[4];
	////////////////
	int sf[2];	//<! 隣接要素のIndex
	int gf[2];
	unsigned int df[2];
};
 */

// ----------------------------------

//! 四面体分割の整合性をチェック
bool CheckTet(
    const std::vector<CDynTet>& tet,
    const std::vector<CDynPointTet>& vertex);

//! 四面体分割の整合性をチェック
bool CheckTet(const std::vector<CDynTet>& tet);

//! ３角形分割の整合性をチェック
//bool CheckTri(const std::vector<STri3D>& tri);

////////////////////////////////////////////////////////////////


//! 四面体分割の辺のリストを得る
bool MakeEdgeTet
 (unsigned int& nedge,
  unsigned int*& edge_ind,
  unsigned int*& edge,
  const std::vector<CDynTet>& tet,
  const unsigned int nnode);

//! 各ノードを囲む四面体の１つを作る
bool MakeOneTetSurNo
 (const std::vector<CDynTet>& tet,
  std::vector<CDynPointTet>& point);



//! 四面体を囲む四面体を作る
//bool MakeTetSurTet(std::vector<CDynTet>& tet);


//! ６面体を囲む６面体を作る
//bool MakeHexSurHex(std::vector<SHex>& aHex);
//! ３角形を囲む３角形を作る
//bool MakeTriSurTri(std::vector<STri3D>& tri);

////////////////
/*
//! 点を囲む６面体のリストを作る
bool MakeHexSurNo(unsigned int& ntetsuno,
				  unsigned int*& tetsuno_ind,
				  unsigned int*& tetsuno,
				  const std::vector<SHex>& aHex,
				  const unsigned int nnode );
 */

/*
//! 点を囲む四面体のリストを作る
bool MakeTetSurNo
 (const std::vector<CDynTet>& aTet,
  const unsigned int nnode,
  unsigned int& ntetsuno,
  unsigned int*& tetsuno_ind,
  unsigned int*& tetsuno );
  */

/*
//! 点を囲む３角形のリストを作る
bool MakeTriSurNo(unsigned int& ntrisuno,
				  unsigned int*& trisuno_ind,
				  unsigned int*& trisuno,
				  const std::vector<STri3D>& tri,
				  const unsigned int nnode );
 */

////////////////

/*
//! 四面体を囲む３角形を作る
bool MakeTriSurTet(std::vector<STet>& in_tet,
				   std::vector<STri3D>& in_tri,
				   const std::vector<CPoint3D>& vertex );
 */

////////////////
/*
//! 四面体分割の境界の３角形を作る
bool MakeOuterBoundTet( const std::vector<STet  >& aTet, std::vector<STri3D>& aTri );
 */
//! ６面体分割の境界の４角形を作る
//bool MakeOuterBoundTri( const std::vector<STri3D>& aTri, std::vector<SBar  >& aBar );

////////////////

/*
// 隣合う三角形の単位法線同士がcrt_dot以下ならそれは切り離す
bool MakeColorCoding_Tri3D( const std::vector<STri3D>& aTri, const std::vector<CVector3>& aVec3D, 
						   std::vector<int>& aColorTri, unsigned int& nColor, double dot_crt);
 */

/*
bool MakeBar_fromColorCodingTri( const std::vector<STri3D>& aTri, 
								const std::vector<int>& aColorTri, unsigned int nColor, std::vector<SBar>& aBar);
 */

/*
bool MakeColorCodingBar( const std::vector<SBar> aBar, 
							 std::vector<int>& aColorBar, unsigned int& nColorBar);

bool MakeColorCodingBar( const std::vector<SBar>& aBar,  const std::vector<CVector3>& aVec, 
						  std::vector<int>& aColorBar, unsigned int& nColor);
*/						  


////////////////////////////////////////////////////////////////


/*
//! ３角形の辺に点を加える
bool AddPointTri_Edge(const unsigned int ipo_ins,
					  const unsigned int itr0,
					  const unsigned int iedtri0,
					  std::vector<STri3D>& aTri );
*/
/*
//! ３角形の面に点を加える
bool AddPointTri_Face(const unsigned int ipo_ins,
					  const unsigned int itr0,
					  std::vector<STri3D>& aTri );q
 */
/*
bool GetAddPointEdgeCrt(const ElemAroundEdge& elared,
						const CVector3& add_point,
						double& max_crt,
						const std::vector<CPoint3D>& node);
 */

/*
//! 四面体の辺に点を加える
bool AddPointTet_Edge
 (const ElemAroundEdge& elared,
  const unsigned int ipo_ins,
  std::vector<CEPo3D>& node,
  std::vector<CDynTet>& tet);
 */

/**
 * @breaf Add point inside tetrahedra and maintain delaunay
 * @param tmp_buffer should be an array of -1. (if input values are all -1, then output values are -1)
 */
void AddPointTetDelaunay(
  unsigned int ip_ins,
  unsigned int itet_ins,
  std::vector<CDynPointTet>& aPo3D,
  std::vector<CDynTet>& aSTet,
  std::vector<CVec3d>& aCent,
  std::vector<int>& tmp_buffer);


/*
//! 四面体の中に点を加える
bool AddPointTet_Elem
 (const unsigned int itet0,
  const unsigned int ipo_ins,
  std::vector<CEPo3D>& po,
  std::vector<CDynTet>& tet);
 */

//! 四面体の面に点を加える
/*
bool AddPointTet_Face
 (const unsigned int itet0,
  const unsigned int ifatet,
  const unsigned int ipo_ins,
  std::vector<CEPo3D>& po,
  std::vector<CDynTet>& tet);
 */

/*
//! 辺を消去する
bool EdgeDel(const ElemAroundEdge& elared,
             const ElemAroundPoint& elarpo0,
             std::vector<STet>& tet,
             std::vector<CPoint3D>& node);
*/
/*
bool GetEdgeDelCrt(const ElemAroundEdge& elared,
				   const ElemAroundPoint& elarpo,
				   double& max_crt,
				   const std::vector<STet>& tet,
				   const std::vector<CPoint3D>& node );
 */

//! 辺周りの要素を取得
bool MakeElemAroundEdge
 ( ElemAroundEdge& elared,
  const int itet0,
  const int idedge0,
  const std::vector<CDynTet>& tet );

/*
double MaxCrtElemAroundPoint(const ElemAroundPoint& elarpo,
						  const std::vector<STet>& tet,
						  const std::vector<CPoint3D>& node);
 */

/**
 * @brief get list of tetrahedra around a point
 */
bool MakeElemAroundPoint
 ( ElemAroundPoint& elarpo,
  const unsigned int itet0,
  const unsigned int inoel0,
  const std::vector<CDynTet>& tet );



/*
double MaxCrtElemAroundEdge(const ElemAroundEdge& elared,
						  const std::vector<STet>& tet,
						  const std::vector<CPoint3D>& node);
 */
/*
bool GetEdgeSwapPtnCrt(const ElemAroundEdge& elared,
				int& ptn,
				double& min_max_crt,
				const std::vector<CPoint3D>& node );
 */

/*
bool EdgeSwapTet
 (const ElemAroundEdge& elared,
  const int ptn,
  std::vector<CDynTet>& tet,
  std::vector<CEPo3D>& point );
 */

// -------------------------

/*
bool GetFaceSwapCrt(const int itet0,
					const int iface0,
					double& max_crt,
					const std::vector<CETet>& tet,
					const std::vector<CEPo3D>& node );
 */

/*
bool FaceSwapTet
 (const unsigned int itet0,
  const unsigned int iface0,
  std::vector<CDynTet>& tet,
  std::vector<CEPo3D>& node);
*/

/*
bool DelaunayAroundPointTet(const unsigned int ipoin,
						 std::vector<CEPo3D>& aPo3D,
						 std::vector<CETet>& aTet);
 */

/*
bool DelaunayAroundPointTri(unsigned int itri, 
							unsigned int inotri, 
							std::vector<CVector3>& aVec,
							std::vector<STri3D>& aTri);
*/
/*
bool FlipEdgeTri( const unsigned int itri0, const unsigned int ied0,
			   std::vector<STri3D>& tri );
*/
// ------------------------------

/*
bool ReconnectCap
 (const unsigned int itet0,
  std::vector<CETet>& tet,
  std::vector<CEPo3D>& node);
*/
/*
bool ReconnectSliver
 (const unsigned int itet0,
  std::vector<CETet>& tet,
  std::vector<CEPo3D>& node);
*/
/*
bool Reconnect
 (std::vector<CETet>& tet,
  std::vector<CEPo3D>& node);
*/




// ---------------------------------------------

//! volume of tetrahedra
inline double TetVolume(
    unsigned int iv1,
    unsigned int iv2,
    unsigned int iv3,
    unsigned int iv4,
    const std::vector<CDynPointTet>& point)
{
	return Volume_Tet( point[iv1].p, point[iv2].p, point[iv3].p, point[iv4].p );
}

//! volume of tetrahedra
inline double TetVolume
(const CDynTet& tet,
 const std::vector<CDynPointTet>& node)
{
	return Volume_Tet(
		node[ tet.v[0] ].p,
		node[ tet.v[1] ].p, 
		node[ tet.v[2] ].p, 
		node[ tet.v[3] ].p );
}

//! volume of tetrahedra
inline double TetVolume
(int ielem,
 const std::vector<CDynTet>& tet,
 const std::vector<CDynPointTet>& node)
{
	return TetVolume( tet[ielem], node );
}

////////////////////////////////////////////////

//! ３角形の面積
inline double TriArea(const CDynTet& tet, 
		const int iface, 
		const std::vector<CDynPointTet>& node )
{
	return Area_Tri(
        node[ tet.v[(int)noelTetFace[iface][0]] ].p,
        node[ tet.v[(int)noelTetFace[iface][1]] ].p,
        node[ tet.v[(int)noelTetFace[iface][2]] ].p );
}

//! ３角形の面積
inline double TriArea(const int itet,
		const int iface,
		const std::vector<CDynTet>& tet,
		const std::vector<CDynPointTet>& node )
{
	return TriArea(tet[itet],iface,node);
}

// ----------------

inline double SquareLongestEdgeLength(const int itet,
		const std::vector<CDynPointTet>& node,
		const std::vector<CDynTet>& tet )
{
	return SqareLongestEdgeLength(
		node[tet[itet].v[0]].p,
		node[tet[itet].v[1]].p,
		node[tet[itet].v[2]].p,
		node[tet[itet].v[3]].p );
}

// -------------------------

inline double LongestEdgeLength(const int itet,
		const std::vector<CDynPointTet>& node,
		const std::vector<CDynTet>& tet )
{
	return sqrt( SqareLongestEdgeLength(
			node[tet[itet].v[0]].p,
			node[tet[itet].v[1]].p,
			node[tet[itet].v[2]].p,
			node[tet[itet].v[3]].p ) );
}

// ------------------------

inline double SquareShortestEdgeLength(const int itet,
		const std::vector<CDynPointTet>& node,
		const std::vector<CDynTet>& tet )
{
	return SqareShortestEdgeLength(
		node[tet[itet].v[0]].p,
		node[tet[itet].v[1]].p,
		node[tet[itet].v[2]].p,
		node[tet[itet].v[3]].p );
}

// ------------------------

inline double ShortestEdgeLength(const int itet,
		const std::vector<CDynPointTet>& node,
		const std::vector<CDynTet>& tet )
{
	return sqrt( SqareShortestEdgeLength(
		node[tet[itet].v[0]].p,
		node[tet[itet].v[1]].p,
		node[tet[itet].v[2]].p,
		node[tet[itet].v[3]].p ) );
}

// ------------------------

//! 法線の取得(長さ１とは限らない)
inline void Normal(
		delfem2::CVec3d& vnorm,
		const unsigned int itet0, const unsigned int iface0, 
		const std::vector<CDynTet>& tet, const std::vector<CDynPointTet>& node)
{
	assert( itet0 < tet.size() );
	assert( iface0 < 4 );
	Normal(vnorm,
         node[ tet[itet0].v[(int)noelTetFace[iface0][0] ] ].p,
         node[ tet[itet0].v[(int)noelTetFace[iface0][1] ] ].p,
         node[ tet[itet0].v[(int)noelTetFace[iface0][2] ] ].p );
}

// ------------------------

//! ３角形の単位法線
inline void UnitNormal(
		delfem2::CVec3d& vnorm,
		const int itet0,
		const int iface0,
		const std::vector<CDynTet>& aTet,
		const std::vector<CDynPointTet>& aPo )
{
	UnitNormal(vnorm,
             aPo[ aTet[itet0].v[ (int)noelTetFace[iface0][0] ] ].p,
             aPo[ aTet[itet0].v[ (int)noelTetFace[iface0][1] ] ].p,
             aPo[ aTet[itet0].v[ (int)noelTetFace[iface0][2] ] ].p );
}

/*
//! ３角形の単位法線
inline void UnitNormal(
		CVector3& vnorm,
		const int itri,
		const std::vector<STri3D>& aTri,
		const std::vector<CVector3>& aVec )
{
  UnitNormal(vnorm,
             aVec[ aTri[itri].v[0] ],
             aVec[ aTri[itri].v[1] ],
             aVec[ aTri[itri].v[2] ] );
}
 */

////////////////////////////////////////////////

//! ３角形の外周円
inline double Circumradius(
		const int itet0, 
		const std::vector<CDynPointTet>& node,
		const std::vector<CDynTet>& tet )
{
	return sqrt( SquareCircumradius(
		node[ tet[itet0].v[0] ].p,
		node[ tet[itet0].v[1] ].p,
		node[ tet[itet0].v[2] ].p,
		node[ tet[itet0].v[3] ].p) );
}



////////////////////////////////////////////////

inline double Criterion_Asp(
		const delfem2::CVec3d& po0,
		const delfem2::CVec3d& po1,
		const delfem2::CVec3d& po2,
		const delfem2::CVec3d& po3)
{
	double saface = Area_Tri( po1, po3, po2 )
				+	Area_Tri( po0, po2, po3 )
				+	Area_Tri( po0, po3, po1 )
				+ Area_Tri( po0, po1, po2 );
	double inscribed_radius = Volume_Tet(po0,po1,po2,po3) * 3.0 / saface;
	double circum_radius = Circumradius(po0,po1,po2,po3);
	return circum_radius / inscribed_radius;
}

inline double Criterion_Asp(const CDynTet& tet,
		const std::vector<CDynPointTet>& node)
{
	return Criterion_Asp(
		node[tet.v[0]].p,
		node[tet.v[1]].p,
		node[tet.v[2]].p,
		node[tet.v[3]].p);
}

inline double Criterion_Asp(int ielem,
		const std::vector<CDynPointTet>& node,
		const std::vector<CDynTet>& tet)
{
	return Criterion_Asp(tet[ielem],node);
}

////////////////////////////////////////////////
/*
inline double Criterion_PLJ(
		const CVector3& po0, 
		const CVector3& po1, 
		const CVector3& po2, 
		const CVector3& po3)
{
	const double longest_edge = LongestEdgeLength(po0,po1,po2,po3);
	double saface = TriArea( po1, po3, po2 )
				+	TriArea( po0, po2, po3 )
				+	TriArea( po0, po3, po1 )
				+   TriArea( po0, po1, po2 );
	return (saface*longest_edge)/(TetVolume(po0,po1,po2,po3)*3.0);
}

inline double Criterion_PLJ(
		const STet& tet,
		const std::vector<CPoint3D>& node)
{
	return Criterion_PLJ(
		node[tet.v[0]].p,
		node[tet.v[1]].p,
		node[tet.v[2]].p,
		node[tet.v[3]].p );
}

inline double Criterion_PLJ(
		int ielem,
		const std::vector<CPoint3D>& node,
		const std::vector<STet>& tet)
{
	return Criterion_PLJ(tet[ielem],node);
}

inline double Criterion_PLJ(
		const CPoint3D& po0, 
		const CPoint3D& po1, 
		const CPoint3D& po2, 
		const CPoint3D& po3)
{
	return Criterion_PLJ(po0.p, po1.p, po2.p, po3.p);
}
 */

////////////////////////////////////////////////

/*! 
@brief ドロネー条件のチェック
v1, v2, v3の外接球の中でもっとも半径の大きな球(中心はv1,v2,v3の外心）の中にv4が含まれているかどうか調べる
@param[in] v1 点
@param[in] v2 点
@param[in] v3 点
@param[in] v4 点
@retval 1 入っていなければ
@retval 0 微妙なら
@retval -1 入っていれば
*/
inline int DetDelaunay
(const delfem2::CVec3d& v1,
 const delfem2::CVec3d& v2,
 const delfem2::CVec3d& v3,
 const delfem2::CVec3d& v4)
{
	// ３角形v1,v2,v3の外接円の中心を求める。
	const double dtmp1 = (v2.x-v3.x)*(v2.x-v3.x)+(v2.y-v3.y)*(v2.y-v3.y)+(v2.z-v3.z)*(v2.z-v3.z);
	const double dtmp2 = (v3.x-v1.x)*(v3.x-v1.x)+(v3.y-v1.y)*(v3.y-v1.y)+(v3.z-v1.z)*(v3.z-v1.z);
	const double dtmp3 = (v1.x-v2.x)*(v1.x-v2.x)+(v1.y-v2.y)*(v1.y-v2.y)+(v1.z-v2.z)*(v1.z-v2.z);

	double qarea = SquareTriArea(v1,v2,v3);
	const double etmp1 = dtmp1*(dtmp2+dtmp3-dtmp1) / (16.0 * qarea );
	const double etmp2 = dtmp2*(dtmp3+dtmp1-dtmp2) / (16.0 * qarea );
	const double etmp3 = dtmp3*(dtmp1+dtmp2-dtmp3) / (16.0 * qarea );

	delfem2::CVec3d out_center(
		etmp1*v1.x + etmp2*v2.x + etmp3*v3.x,
		etmp1*v1.y + etmp2*v2.y + etmp3*v3.y,
		etmp1*v1.z + etmp2*v2.z + etmp3*v3.z );
	
	const double qradius = SquareDistance( out_center, v1 );
	assert( fabs( qradius - SquareDistance(out_center,v2) ) < 1.0e-10 );
	assert( fabs( qradius - SquareDistance(out_center,v3) ) < 1.0e-10 );
	assert( fabs( Height(v1,v2,v3,out_center) ) < 1.0e-10 );

	const double distance = SquareDistance( out_center, v4 );
	if( distance > qradius + 1.0e-6 ){
		return 1;
	}
	else{
		if( distance < qradius - 1.0e-6 ){
			return -1;
		}
		else{
			return 0;
		}
	}
	assert(0);
	return -1;
}

/*! 
@brief ドロネー条件のチェック
v1, v2, v3，v4の外接球の中にv5が含まれているかどうか調べる
@param[in] v1 点
@param[in] v2 点
@param[in] v3 点
@param[in] v4 点
@param[in] v5 点
@retval >0 入っていなければ
@retval <0 入っていれば
*/
inline double DetDelaunay3D
(const delfem2::CVec3d& v1,
 const delfem2::CVec3d& v2,
 const delfem2::CVec3d& v3,
 const delfem2::CVec3d& v4,
 const delfem2::CVec3d& v5)
{
	const double a[12] = {
		v1.x-v5.x,	//  0
		v2.x-v5.x,	//  1
		v3.x-v5.x,	//  2
		v4.x-v5.x,	//  3
		v1.y-v5.y,	//  4
		v2.y-v5.y,	//  5
		v3.y-v5.y,	//  6
		v4.y-v5.y,	//  7
		v1.z-v5.z,	//  8
		v2.z-v5.z,	//  9
		v3.z-v5.z,	// 10
		v4.z-v5.z,	// 11
	};
	const double b[6] = {
		a[ 6]*a[11]-a[ 7]*a[10],	// 0
		a[ 5]*a[11]-a[ 7]*a[ 9],	// 1
		a[ 5]*a[10]-a[ 6]*a[ 9],	// 2
		a[ 7]*a[ 8]-a[ 4]*a[11],	// 3
		a[ 6]*a[ 8]-a[ 4]*a[10],	// 4
		a[ 4]*a[ 9]-a[ 5]*a[ 8],	// 5 
	};
	return
      -( a[0]*(v1.x+v5.x)+a[4]*(v1.y+v5.y)+a[ 8]*(v1.z+v5.z) )*( a[ 1]*b[0]-a[ 2]*b[1]+a[ 3]*b[2] )
			+( a[1]*(v2.x+v5.x)+a[5]*(v2.y+v5.y)+a[ 9]*(v2.z+v5.z) )*( a[ 0]*b[0]+a[ 2]*b[3]-a[ 3]*b[4] )
			-( a[2]*(v3.x+v5.x)+a[6]*(v3.y+v5.y)+a[10]*(v3.z+v5.z) )*( a[ 0]*b[1]+a[ 1]*b[3]+a[ 3]*b[5] )
			+( a[3]*(v4.x+v5.x)+a[7]*(v4.y+v5.y)+a[11]*(v4.z+v5.z) )*( a[ 0]*b[2]+a[ 1]*b[4]+a[ 2]*b[5] );
}

/*! 
@brief ドロネー条件のチェック
v1, v2, v3，v4の外接球の中にv5が含まれているかどうか調べる
@param[in] v1 点のIndex
@param[in] v2 点のIndex
@param[in] v3 点のIndex
@param[in] v4 点のIndex
@param[in] v5 点のIndex
@param[in] node 座標配列
@retval >0 入っていなければ
@retval <0 入っていれば
*/
inline double DetDelaunay3D
( const int v1, const int v2, const int v3, const int v4, const int v5,
 const std::vector<delfem2::CVec3d>& node )
{
	return DetDelaunay3D(node[v1],node[v2],node[v3],node[v4],node[v5]);
}

} // namespace delfem2


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/dtet_v3.cpp"
#endif

#endif
