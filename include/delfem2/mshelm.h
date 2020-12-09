/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file functions to analyze mesh topology for mixed meshes
 * @details the functions only care about the topology. Geometry (coordinate) information is not handled in this file
 */

// TODO: separate mixed elem

#ifndef DFM2_MSHELM_H
#define DFM2_MSHELM_H

#include "delfem2/dfm2_inline.h"

namespace delfem2
{

enum MESHELEM_TYPE
{
  MESHELEM_TRI = 0,
  MESHELEM_TET = 1,
  MESHELEM_QUAD = 2,
  MESHELEM_HEX = 3,
  MESHELEM_PYRAMID = 4,
  MESHELEM_WEDGE = 5,
  MESHELEM_VOX = 6,
  MESHELEM_LINE = 7,
};

// 5 : VTK_TRIANGLE
// 9 : VTK_QUAD
// 10: VTK_TETRA
// 12: VTK_HEXAHEDRON
// 13: VTK_WEDGE
// 14: VTK_PYRAMD

////
const int mapMeshElemType2NNodeElem[8] = {
  3, // TRI 0
  4, // TET 1
  4, // QUAD 2
  8, // HEX 3
  5, // PYRAM 4
  6, // WEDGE 5
  8,  // VOX 6
  2, // LINE 7
};
////
const int mapMeshElemType2NFaceElem[7] = {
  3, // TRI
  4, // TET
  4, // QUAD
  6, // HEX
  5, // Pyramid
  5, // Wedge
  6  // VOX
};
////
const int mapMeshElemTypeIndFace2NNoelElemFace[7][6] =
{
  {2,2,2,-1,-1,-1}, // TRI
  {3,3,3, 3,-1,-1}, // TET
  {2,2,2, 2,-1,-1}, // QUAD
  {4,4,4, 4, 4, 4}, // HEX
  {4,3,3, 3, 3,-1}, // Pyramid
  {4,4,4, 3, 3,-1}, // Wedge
  {4,4,4, 4, 4, 4}, // VOX
};
////
const int noelElemFace_Tri[3][4] = {
  { 1, 2,-1,-1 }, //edge 0
  { 2, 0,-1,-1 }, //edge 1
  { 0, 1,-1,-1 }, //edge 2
};
const int noelElemFace_Tet[4][4] ={
  { 1, 2, 3,-1 },
  { 0, 3, 2,-1 },
  { 0, 1, 3,-1 },
  { 0, 2, 1,-1 }
};
const int noelElemFace_Quad[4][4] = {
  { 0, 1,-1,-1 },
  { 1, 2,-1,-1 },
  { 2, 3,-1,-1 },
  { 3, 0,-1,-1 },
};
const int noelElemFace_Vox[6][4] = { // this numbering is corresponds to VTK_VOX
  { 0, 4, 6, 2 }, // -x
  { 1, 3, 7, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 2, 6, 7, 3 }, // +y
  { 0, 2, 3, 1 }, // -z
  { 4, 5, 7, 6 }  // +z
};
const int noelElemFace_Hex[6][4] = { // this numbering is corresponds to VTK_HEX
  { 0, 4, 7, 3 }, // -x
  { 1, 2, 6, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 3, 7, 6, 2 }, // +y
  { 0, 3, 2, 1 }, // -z
  { 4, 5, 6, 7 }  // +z
};
const int noelElemFace_Pyram[5][4] = {
  { 0, 3, 2, 1},
  { 0, 1, 4,-1},
  { 1, 2, 4,-1},
  { 2, 3, 4,-1},
  { 3, 0, 4,-1}
};
const int noelElemFace_Wed[5][4] = {
  { 0, 2, 5, 3},
  { 1, 0, 3, 4},
  { 2, 1, 4, 5},
  { 0, 1, 2,-1},
  { 3, 5, 4,-1}
};

const int mapMeshElemType2NEdgeElem[7] = {
  3, // TRI
  6, // TET
  4, // QUAD
  12, // HEX
  8, // PYRAM
  9, // WEDGE
  12  // VOX
};
const int noelElemEdge_Tri[3][2] = {
  { 0, 1},
  { 1, 2},
  { 2, 0},
};
const int noelElemEdge_Tet[6][2] = {
  { 0, 1},
  { 0, 2},
  { 0, 3},
  { 1, 2},
  { 1, 3},
  { 2, 2},
};
const int noelElemEdge_Quad[4][2] = {
  { 0, 1},
  { 1, 2},
  { 2, 3},
  { 3, 0},
};
const int noelElemEdge_Hex[12][2] = {
  {0,1},{1,2},{2,3},{3,0},
  {4,5},{5,6},{6,7},{7,4},
  {0,4},{1,5},{2,6},{3,7}
};
// TODO: this looks wrong...
const int noelElemEdge_Vox[12][2] = {
  {0,1},{3,2},{4,5},{7,6},
  {0,3},{1,2},{4,7},{5,6},
  {0,4},{1,5},{3,7},{2,6}
};

// -----------------------------
inline int nNodeElem(MESHELEM_TYPE type){
  return mapMeshElemType2NNodeElem[type];
}
inline int nFaceElem(MESHELEM_TYPE type){
  return mapMeshElemType2NFaceElem[type];
}
inline int nNodeElemFace(MESHELEM_TYPE type,int iface){
  return mapMeshElemTypeIndFace2NNoelElemFace[type][iface];
}
inline const int (*noelElemFace(MESHELEM_TYPE type))[4]
{
  const int (*noelElemFace0[7])[4] = {
    noelElemFace_Tri,
    noelElemFace_Tet,
    noelElemFace_Quad,
    noelElemFace_Hex,
    noelElemFace_Pyram,
    noelElemFace_Wed,
    noelElemFace_Vox,
  };
  return noelElemFace0[type];
}
inline const int (*noelElemEdge(MESHELEM_TYPE type))[2]
{
  const int (*noelElemEdge0[4])[2] = {
    noelElemEdge_Tri,
    noelElemEdge_Tet,
    noelElemEdge_Quad,
    noelElemEdge_Hex
  };
  return noelElemEdge0[type];
}


}


#endif
