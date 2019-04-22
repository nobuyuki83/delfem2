#ifndef DYNTRI_H
#define DYNTRI_H

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

class CEPo2{
public:
	CEPo2(){
    e = -1;
    d = 0;
  }
	CEPo2( const CEPo2& rhs )
		: e(rhs.e), d(rhs.d), p(rhs.p) {}
  CEPo2(const CVector3 p, int ielem, unsigned int idir)
    : e(ielem), d(idir), p(p){}
  CEPo2(const CVector3& p) :
    e(-1), d(0), p(p){}
public:
  int e;  
  int d; 
  ////
  CVector3 p; // position
};

bool MakePointSurTri(const std::vector<ETri>& aTri, const unsigned int npoin,
                     unsigned int* const elsup_ind, unsigned int& nelsup, unsigned int*& elsup );

void MakeEdge(unsigned int* const edge_ind, unsigned int& nedge, unsigned int*& edge,
              const std::vector<ETri>& aTri, const unsigned int npoin,
              const unsigned int* elsup_ind, const unsigned int nelsup, const unsigned int* elsup);

bool MakeInnerRelationTri(std::vector<ETri>& aTri, const unsigned int npoin,
                          const unsigned int* elsup_ind, const unsigned int nelsup, const unsigned int* elsup);

bool CheckTri(const std::vector<ETri>& aTri);


bool FindEdge(int& itri0, int& inotri0, int& inotri1,
              ///
              const int ipo0, const int ipo1,
              const std::vector<CEPo2>& aPo,
              const std::vector<ETri>& aTri);


///////////////////////////////////////////////////////////////
// topology edit

bool InsertPoint_ElemEdge(const int ipo_ins,    //the index of the new point
                          const int itri_ins,  //triangle index
                          const int ied_ins,  //edge index
                          std::vector<CEPo2>& aPo,
                          std::vector<ETri>& aTri );

bool InsertPoint_Elem(const int ipo_ins,
                      const int itri_ins,
                      std::vector<CEPo2>& aPo,
                      std::vector<ETri>& aTri);

bool DelaunayAroundPoint(int ipo0,
                         std::vector<CEPo2>& aPo,
                         std::vector<ETri>& aTri);

bool FlipEdge(int itri0, int ied0,
              std::vector<CEPo2>& aPo,
              std::vector<ETri>& aTri);

void MoveCCW(int& itri_cur,
             int& inotri_cur,
             bool& flag_is_wall,
             ////
             std::vector<CEPo2>& aPo,
             std::vector<ETri>& aTri);

bool DelaunayAroundPoint(int ipo0,
                         std::vector<CEPo2>& aPo,
                         std::vector<ETri>& aTri);

bool DeleteTri(unsigned int itri_to,
               std::vector<CEPo2>& aPo,
               std::vector<ETri>& aTri);

bool Collapse_ElemEdge(const int itri_del,
                       const int ied_del,
                       std::vector<CEPo2>& aPo,
                       std::vector<ETri>& aTri);

void GetTriAryAroundPoint(int ipoin,
                          std::vector< std::pair<unsigned int,unsigned int> >& aTriSurPo,
                          const std::vector<CEPo2>& aPo,
                          const std::vector<ETri>& aTri);


/////////////////////////////////////////////////////////////////////////////////////////////////

void extractHoles(std::vector< std::vector<int> >& aIndP_Hole,
                  const int npo,
                  const std::vector<ETri>& aETri);


#endif // #endif SURFACE_MESH_H
