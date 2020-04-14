/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_DTRI_H
#define DFM2_DTRI_H

#include <vector>
#include "delfem2/dfm2_inline.h"

// ------------------------------------------------------------------

namespace delfem2 {

// (6-i-j)%3;
static const unsigned int relTriTri[3][3] = {
	{ 0, 2, 1 }, //  0
	{ 2, 1, 0 }, //  1 
	{ 1, 0, 2 }, //  2
};
static const unsigned int invRelTriTri[3] = { 0, 1, 2 };

class CDynTri{
public:
	int v[3];	//!< index of vertex
	int s2[3];	//!< index of face that is adjacent to ith edge; The 0th edge is the edge facing 0th vertex
  int r2[3];	//!< relationship of vertex index between two adjacent faces
};

class CDynPntSur{
public:
	CDynPntSur(){ e = -1; d = 0; }
	CDynPntSur( const CDynPntSur& rhs )
		: e(rhs.e), d(rhs.d) {}
  CDynPntSur(int ielem, unsigned int idir)
    : e(ielem), d(idir) {}
public:
  int e;  //<! index of elements this can be zero
  unsigned int d;
};

DFM2_INLINE bool JArray_MakeElSuP
 (std::vector<int>& elsup_ind, std::vector<int>& elsup,
  const std::vector<CDynTri>& aTri, const unsigned int npoin);

DFM2_INLINE void JArray_PSuP
 (std::vector<int>& psup_ind, std::vector<int>& psup,
  const std::vector<CDynTri>& aTri, const unsigned int npoin,
  const std::vector<int>& elsup_ind, const std::vector<int>& elsup);

DFM2_INLINE void MakeInnerRelationTri
 (std::vector<CDynTri>& aTri, const unsigned int npoin,
  const std::vector<int>& elsup_ind,
  const std::vector<int>& elsup);

DFM2_INLINE void JArray_PSuP
 (unsigned int* const edge_ind, unsigned int& nedge, unsigned int*& edge,
  const std::vector<CDynTri>& aTri, const unsigned int npoin,
  const std::vector<int>& elsup_ind, const std::vector<int>& elsup);

DFM2_INLINE bool CheckTri
 (const std::vector<CDynTri>& aETri);

DFM2_INLINE bool CheckTri
 (const std::vector<CDynPntSur>& aEPo2,
  const std::vector<CDynTri>& aETri,
  bool is_assert=true);

DFM2_INLINE void InitializeMesh
 (std::vector<CDynPntSur>& aEPo2,
  std::vector<CDynTri>& aETri,
  const unsigned int* aTri,  int nTri,
  int nXYZ);

DFM2_INLINE bool FindEdge_LookAroundPoint
 (unsigned int &itri0,
  unsigned int &inotri0,
  unsigned &inotri1,
  //
  const int ipo0, const int ipo1,
  const std::vector<CDynPntSur>& aPo,
  const std::vector<CDynTri>& aTri);

DFM2_INLINE bool FindEdge_LookAllTriangles
 (int& itri0, int& iedtri0,
  //
  const int ipo0, const int ipo1,
  const std::vector<CDynTri>& aETri);

DFM2_INLINE void GetTriArrayAroundPoint
 (std::vector< std::pair<int,int> >& aTriSurPo,
  int ipoin,
  const std::vector<CDynPntSur>& aEPo2,
  const std::vector<CDynTri>& aETri);

DFM2_INLINE void MoveCCW
 (int& itri_cur,
  unsigned int &inotri_cur,
  bool& flag_is_wall,
  //
  std::vector<CDynTri>& aTri);

// ---------------
// topology edit

DFM2_INLINE bool FlipEdge
 (unsigned int itri0, unsigned int ied0,
  std::vector<CDynPntSur>& aPo,
  std::vector<CDynTri>& aTri);

// ----------------------
// insert point

DFM2_INLINE bool InsertPoint_ElemEdge
 (const int ipo_ins,  //!< the index of the new point
  const int itri_ins, //!< triangle index
  const int ied_ins,  //!< edge index
  std::vector<CDynPntSur>& aEPo2,
  std::vector<CDynTri>& aETri );

DFM2_INLINE bool InsertPoint_Elem
 (const int ipo_ins,
  const int itri_ins,
  std::vector<CDynPntSur>& aEPo2,
  std::vector<CDynTri>& aETri);

// -----------------
// delete point

DFM2_INLINE bool DeleteTri
 (int itri_to,
  std::vector<CDynPntSur>& aEPo2,
  std::vector<CDynTri>& aETri);

DFM2_INLINE bool Collapse_ElemEdge
 (const int itri_del,
  const int ied_del,
  std::vector<CDynPntSur>& aEPo2,
  std::vector<CDynTri>& aETri);


// ------------------------------------------

DFM2_INLINE void extractHoles
 (std::vector< std::vector<int> >& aIndP_Hole,
  const int npo,
  const std::vector<CDynTri>& aETri);

}


#ifdef DFM2_HEADER_ONLY
#  include "delfem2/dtri.cpp"
#endif

#endif // #endif DTRI_H
