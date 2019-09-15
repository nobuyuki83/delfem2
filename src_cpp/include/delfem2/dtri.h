/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DYNTRI_H
#define DYNTRI_H

#include <vector>

//////////////////////////////////////////////////////////////////////////////////////////////////

// (6-i-j)%3;
static const unsigned int relTriTri[3][3] = {
	{ 0, 2, 1 }, //  0
	{ 2, 1, 0 }, //  1 
	{ 1, 0, 2 }, //  2
};
static const unsigned int invRelTriTri[3] = { 0, 1, 2 };

class ETri{
public:
	int v[3];	//!< index of vertex
	int s2[3];	//!< index of face that is adjacent to ith edge; The 0th edge is the edge facing 0th vertex
  int r2[3];	//!< relationship of vertex index between two adjacent faces
};

class CEPo2{
public:
	CEPo2(){ e = -1; d = 0; }
	CEPo2( const CEPo2& rhs )
		: e(rhs.e), d(rhs.d) {}
  CEPo2(int ielem, unsigned int idir)
    : e(ielem), d(idir) {}
public:
  int e;  //<! index of 
  int d;
};

bool JArray_MakeElSuP(std::vector<int>& elsup_ind, std::vector<int>& elsup,
                     const std::vector<ETri>& aTri, const unsigned int npoin);

void MakeInnerRelationTri(std::vector<ETri>& aTri, const unsigned int npoin,
                          const std::vector<int>& elsup_ind,
                          const std::vector<int>& elsup);

void JArray_PSuP(unsigned int* const edge_ind, unsigned int& nedge, unsigned int*& edge,
              const std::vector<ETri>& aTri, const unsigned int npoin,
              const std::vector<int>& elsup_ind, const std::vector<int>& elsup);

bool CheckTri(const std::vector<ETri>& aETri);
bool CheckTri(const std::vector<CEPo2>& aEPo2,
              const std::vector<ETri>& aETri,
              bool is_assert=true);

void InitializeMesh(std::vector<CEPo2>& aEPo2,
                    std::vector<ETri>& aETri,
                    const unsigned int* aTri,  int nTri,
                    int nXYZ);

bool FindEdge_LookAroundPoint(int& itri0, int& inotri0, int& inotri1,
              ///
              const int ipo0, const int ipo1,
              const std::vector<CEPo2>& aEPo2,
              const std::vector<ETri>& aETri);

bool FindEdge_LookAllTriangles(int& itri0, int& iedtri0,
                               ///
                               const int ipo0, const int ipo1,
                               const std::vector<ETri>& aETri);

void GetTriArrayAroundPoint(std::vector< std::pair<int,int> >& aTriSurPo,
                          int ipoin,
                          const std::vector<CEPo2>& aEPo2,
                          const std::vector<ETri>& aETri);

void MoveCCW(int& itri_cur,
             int& inotri_cur,
             bool& flag_is_wall,
             ////
             std::vector<ETri>& aETri);

///////////////////////////////////////////////////////////////
// topology edit

bool FlipEdge(int itri0, int ied0,
              std::vector<CEPo2>& aEPo2,
              std::vector<ETri>& aETri);

////////////////
// insert point

bool InsertPoint_ElemEdge(const int ipo_ins,    //the index of the new point
                          const int itri_ins,  //triangle index
                          const int ied_ins,  //edge index
                          std::vector<CEPo2>& aEPo2,
                          std::vector<ETri>& aETri );

bool InsertPoint_Elem(const int ipo_ins,
                      const int itri_ins,
                      std::vector<CEPo2>& aEPo2,
                      std::vector<ETri>& aETri);

////////////////
// delete point

bool DeleteTri(int itri_to,
               std::vector<CEPo2>& aEPo2,
               std::vector<ETri>& aETri);

bool Collapse_ElemEdge(const int itri_del,
                       const int ied_del,
                       std::vector<CEPo2>& aEPo2,
                       std::vector<ETri>& aETri);


/////////////////////////////////////////////////////////////////////////////////////////////////

void extractHoles(std::vector< std::vector<int> >& aIndP_Hole,
                  const int npo,
                  const std::vector<ETri>& aETri);


//////////////////////////////////////////////////////////////////////////////////////////////////



#endif // #endif SURFACE_MESH_H
