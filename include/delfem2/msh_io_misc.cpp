/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/msh_io_misc.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <cstring>
#include <climits>

namespace delfem2::msh_iomisc {

DFM2_INLINE void UnitNormalAreaTri3D_(
    double *n, double &a,
    const double *v1, const double *v2, const double *v3) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  a = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * 0.5;
  const double invlen = 0.5 / a;
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

// probably std::stroi is safer to use but it is only for C++11
DFM2_INLINE int MyStoi_(const std::string &str) {
  char *e;
  long d = std::strtol(str.c_str(), &e, 0);
  return (int) d;
}

DFM2_INLINE double MyStod_(const std::string &str) {
  char *e;
  double d = std::strtod(str.c_str(), &e);
  return d;
}

DFM2_INLINE void ParseVtxObj_(
    int &ip,
    int &it,
    int &in,
    char *str) {
  ip = -1;
  it = -1;
  in = -1;
  std::size_t n = strlen(str);
  int icnt = 0;
  int aloc[3];
  for (unsigned int i = 0; i < n; ++i) {
    if (str[i] != '/') { continue; }
    str[i] = '\0';
    aloc[icnt] = i;
    icnt++;
  }
  ip = std::stoi(str);
//  sscanf(str,"%d",&ip);
  ip--;
  if (icnt == 0) {
    return;
  }
  if (icnt == 1) {
    it = std::stoi(str + aloc[0] + 1);
//    sscanf(str+aloc[0]+1,"%d",&it);
    it--;
    return;
  }
  if (icnt == 2) {
    if (aloc[1] - aloc[0] != 1) {
      it = std::stoi(str + aloc[0] + 1);
//      sscanf(str+aloc[0]+1,"%d",&it);
      it--;
    }
    in = std::stoi(str + aloc[1] + 1);
//    sscanf(str+aloc[1]+1,"%d",&in);
    in--;
    return;
  }
}

}

// ----------------------------------------------------

DFM2_INLINE void delfem2::Write_STL(
    const std::string &str,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri) {
  std::ofstream fout(str.c_str(), std::ofstream::out);
  fout << "solid \"hoge\"" << std::endl;
  for (int itri = 0; itri < (int) aTri.size() / 3; itri++) {
    int i0 = aTri[itri * 3 + 0];
    int i1 = aTri[itri * 3 + 1];
    int i2 = aTri[itri * 3 + 2];
    double p0[3] = {aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]};
    double p1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
    double p2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
    double n[3];
    double a;
    msh_iomisc::UnitNormalAreaTri3D_(n, a, p0, p1, p2);
    fout << "  facet normal " << n[0] << " " << n[1] << " " << n[2] << std::endl;
    fout << "    outer loop" << std::endl;
    fout << "      vertex " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
    fout << "      vertex " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
    fout << "      vertex " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
    fout << "    endloop" << std::endl;
    fout << "  endfacet" << std::endl;
  }
  fout << "endsolid" << std::endl;
}

// ------------------------------------------------------------------------

DFM2_INLINE void delfem2::Read_MeshTet3D_TetGen(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<int> &aTet,
    std::vector<int> &aTri) {
  //    std::string fname("example.3");
  //    std::string fname("bunny.1");
  //    std::string fname("bunny_c.1");
  //  std::string fname("camel_c.1");
  //    std::string fname("camel_10000.1");
  {
    std::ifstream fin((fname + std::string(".node")).c_str());
    unsigned int nnode;
    unsigned int i1, i2, i3;
    fin >> nnode >> i1 >> i2 >> i3;
    aXYZ.resize(nnode * 3);
    for (unsigned int ino = 0; ino < nnode; ino++) {
      fin >> i1 >> aXYZ[ino * 3 + 0] >> aXYZ[ino * 3 + 1] >> aXYZ[ino * 3 + 2];
      //        std::cout << i1 << " " << aXYZ[ino*3+0] << " " << aXYZ[ino*3+1] << " " << aXYZ[ino*3+2] << std::endl;
    }
  }
  {
    std::ifstream fin((fname + std::string(".ele")).c_str());
    unsigned int ntet;
    unsigned int i1, i2;
    fin >> ntet >> i1 >> i2;
    aTet.resize(ntet * 4);
    for (unsigned int itet = 0; itet < ntet; itet++) {
      fin >> i1 >> aTet[itet * 4 + 0] >> aTet[itet * 4 + 1] >> aTet[itet * 4 + 2] >> aTet[itet * 4 + 3];
      //        std::cout << i1 << " " << aTet[itet*4+0] << " ";
      //        std::cout << aTet[itet*4+1] << " ";
      //        std::cout << aTet[itet*4+2] << " ";
      //        std::cout << aTet[itet*4+3] << std::endl;
      //        aTet[itet*4+0] -= 1;
      //        aTet[itet*4+1] -= 1;
      //        aTet[itet*4+2] -= 1;
      //        aTet[itet*4+3] -= 1;

    }
  }
  {
    std::ifstream fin((fname + std::string(".face")).c_str());
    int ntri;
    unsigned int i1;
    fin >> ntri >> i1;
    aTri.resize(ntri * 3);
    //      char buff[256];
    for (int itri = 0; itri < ntri; itri++) {
      /*
       int i2;
       fin.getline(buff, 256);
       std::stringstream ss(buff);
       ss>>i1>>aTri[itri*3+0]>>aTri[itri*3+2]>>aTri[itri*3+1]>>i2;
       //        fin >> i1 >> aTri[itri*3+0] >> aTri[itri*3+2] >> aTri[itri*3+1] >> i2;
       */
      fin >> i1 >> aTri[itri * 3 + 0] >> aTri[itri * 3 + 2] >> aTri[itri * 3 + 1];
      //        aTri[itri*3+0] -= 1;
      //        aTri[itri*3+1] -= 1;
      //        aTri[itri*3+2] -= 1;

    }
  }
}

// --------------------------------------------

DFM2_INLINE void delfem2::Write_MeshTri3D(
    std::ofstream &fout,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri,
    int ndim) {
  const int nXYZ = (int) aXYZ.size() / ndim;
  fout << nXYZ << " " << aTri.size() / 3 << std::endl;
  for (int ixyz = 0; ixyz < nXYZ; ++ixyz) {
    for (int idim = 0; idim < ndim; ++idim) {
      fout << aXYZ[ixyz * ndim + idim] << " ";
    }
    fout << std::endl;
  }
  for (int itri = 0; itri < (int) aTri.size() / 3; ++itri) {
    fout << aTri[itri * 3 + 0] << " " << aTri[itri * 3 + 1] << " " << aTri[itri * 3 + 2] << std::endl;
  }
}

DFM2_INLINE void delfem2::Read_MeshTri3D(
    std::ifstream &fin,
    std::vector<double> &aXYZ,
    std::vector<int> &aTri) {
  int nxyz, ntri;
  fin >> nxyz >> ntri;
  aXYZ.resize(nxyz * 3);
  aTri.resize(ntri * 3);
  for (int ixyz = 0; ixyz < nxyz; ++ixyz) {
    fin >> aXYZ[ixyz * 3 + 0] >> aXYZ[ixyz * 3 + 1] >> aXYZ[ixyz * 3 + 2];
  }
  for (int itri = 0; itri < ntri; ++itri) {
    fin >> aTri[itri * 3 + 0] >> aTri[itri * 3 + 1] >> aTri[itri * 3 + 2];
  }
}

DFM2_INLINE void delfem2::WriteVTK_Points(
    std::ofstream &fout,
    const std::string &name,
    const double *aXYZ,
    int nXYZ,
    int nDim) {
  fout << "# vtk DataFile Version 2.0" << std::endl;
  fout << name << std::endl;
  fout << "ASCII" << std::endl;
  fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
  const int np = nXYZ;
  fout << "POINTS " << np << " float" << std::endl;
  if (nDim == 3) {
    for (int ip = 0; ip < np; ++ip) {
      fout << aXYZ[ip * 3 + 0] << " " << aXYZ[ip * 3 + 1] << " " << aXYZ[ip * 3 + 2] << std::endl;
    }
  } else if (nDim == 2) {
    for (int ip = 0; ip < np; ++ip) {
      fout << aXYZ[ip * 2 + 0] << " " << aXYZ[ip * 2 + 1] << " " << 0.0 << std::endl;
    }
  } else {
    assert(0);
  }
}

// 5 : VTK_TRIANGLE
// 9 : VTK_QUAD
// 10: VTK_TETRA
// 12: VTK_HEXAHEDRON
// 13: VTK_WEDGE
// 14: VTK_PYRAMD
void delfem2::WriteVTK_Cells(
    std::ofstream &fout,
    int vtk_elem_type,
    const int *aElem,
    const int nElem) {
  ////
  int nnoel = -1;
  if (vtk_elem_type == 12) { nnoel = 8; }
  else if (vtk_elem_type == 5) { nnoel = 3; }
  else if (vtk_elem_type == 10) { nnoel = 4; }
  else if (vtk_elem_type == 9) { nnoel = 4; }
  const int nelem = nElem;
  fout << "CELLS " << nelem << " " << nelem * (nnoel + 1) << std::endl;
  if (vtk_elem_type == 12) {
    for (int ie = 0; ie < nelem; ++ie) {
      fout << nnoel << " ";
      fout << aElem[ie * 8 + 0] << " " << aElem[ie * 8 + 1] << " " << aElem[ie * 8 + 2] << " " << aElem[ie * 8 + 3]
           << " ";
      fout << aElem[ie * 8 + 4] << " " << aElem[ie * 8 + 5] << " " << aElem[ie * 8 + 6] << " " << aElem[ie * 8 + 7]
           << std::endl;
    }
    fout << "CELL_TYPES " << nelem << std::endl;
    for (int ie = 0; ie < nelem; ++ie) { fout << "12" << std::endl; }
  } else if (vtk_elem_type == 10) {
    for (int ie = 0; ie < nelem; ++ie) {
      fout << nnoel << " ";
      fout << aElem[ie * 4 + 0] << " " << aElem[ie * 4 + 1] << " " << aElem[ie * 4 + 2] << " " << aElem[ie * 4 + 3]
           << std::endl;
    }
    fout << "CELL_TYPES " << nelem << std::endl;
    for (int ie = 0; ie < nelem; ++ie) { fout << "10" << std::endl; }
  } else if (vtk_elem_type == 5) {
    for (int ie = 0; ie < nelem; ++ie) {
      fout << nnoel << " ";
      fout << aElem[ie * 3 + 0] << " " << aElem[ie * 3 + 1] << " " << aElem[ie * 3 + 2] << std::endl;
    }
    fout << "CELL_TYPES " << nelem << std::endl;
    for (int ie = 0; ie < nelem; ++ie) { fout << "5" << std::endl; }
  } else if (vtk_elem_type == 9) {
    for (int ie = 0; ie < nelem; ++ie) {
      fout << nnoel << " ";
      fout << aElem[ie * 4 + 0] << " " << aElem[ie * 4 + 1] << " " << aElem[ie * 4 + 2] << " " << aElem[ie * 4 + 3]
           << std::endl;
    }
    fout << "CELL_TYPES " << nelem << std::endl;
    for (int ie = 0; ie < nelem; ++ie) { fout << "9" << std::endl; }
  }
}

DFM2_INLINE void delfem2::WriteVTK_Cells(
    std::ofstream &fout,
    const std::vector<int> &aTet,
    const std::vector<int> &aPyrm,
    const std::vector<int> &aPrsm) {
  const std::size_t ntet = aTet.size() / 4;
  const std::size_t npyrm = aPyrm.size() / 5;
  const std::size_t nprsm = aPrsm.size() / 6;
  std::size_t nelem = ntet + npyrm + nprsm;
  fout << "CELLS " << nelem << " " << ntet * 5 + npyrm * 6 + nprsm * 7 << std::endl;
  for (std::size_t it = 0; it < ntet; ++it) {
    fout << 4 << " ";
    fout << aTet[it * 4 + 0] << " ";
    fout << aTet[it * 4 + 1] << " ";
    fout << aTet[it * 4 + 2] << " ";
    fout << aTet[it * 4 + 3] << std::endl;
  }
  for (std::size_t ipyrm = 0; ipyrm < npyrm; ++ipyrm) {
    fout << 5 << " ";
    fout << aPyrm[ipyrm * 5 + 0] << " ";
    fout << aPyrm[ipyrm * 5 + 1] << " ";
    fout << aPyrm[ipyrm * 5 + 2] << " ";
    fout << aPyrm[ipyrm * 5 + 3] << " ";
    fout << aPyrm[ipyrm * 5 + 4] << std::endl;
  }
  for (std::size_t iprsm = 0; iprsm < nprsm; ++iprsm) {
    fout << 6 << " ";
    fout << aPrsm[iprsm * 6 + 0] << " ";
    fout << aPrsm[iprsm * 6 + 1] << " ";
    fout << aPrsm[iprsm * 6 + 2] << " ";
    fout << aPrsm[iprsm * 6 + 3] << " ";
    fout << aPrsm[iprsm * 6 + 4] << " ";
    fout << aPrsm[iprsm * 6 + 5] << std::endl;
  }
  fout << "CELL_TYPES " << nelem << std::endl;
  for (std::size_t ie = 0; ie < ntet; ++ie) { fout << "10" << std::endl; }
  for (std::size_t ie = 0; ie < npyrm; ++ie) { fout << "14" << std::endl; }
  for (std::size_t ie = 0; ie < nprsm; ++ie) { fout << "13" << std::endl; }
}

DFM2_INLINE void delfem2::WriteVTK_Data_PointVec(
    std::ofstream &fout,
    const double *aVal,
    int np,
    int nStrideVal,
    int ndim) {
  fout << "VECTORS point_vectors float" << std::endl;
  if (ndim == 3) {
    for (int ip = 0; ip < np; ++ip) {
      fout << aVal[ip * nStrideVal + 0] << " ";
      fout << aVal[ip * nStrideVal + 1] << " ";
      fout << aVal[ip * nStrideVal + 2] << std::endl;
    }
  }
  if (ndim == 2) {
    for (int ip = 0; ip < np; ++ip) {
      fout << aVal[ip * nStrideVal + 0] << " ";
      fout << aVal[ip * nStrideVal + 1] << " ";
      fout << 0.0 << std::endl;
    }
  }
}

DFM2_INLINE void delfem2::WriteVTK_Data_PointScalar(
    std::ofstream &fout,
    const double *aVal,
    int np,
    int nStrideVal) {
  fout << "SCALARS point_scalars float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for (int ip = 0; ip < np; ++ip) {
    fout << aVal[ip * nStrideVal] << std::endl;
  }
}

/*
 void WriteVTK_HexVec
 (const std::string& fpath,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aHex,
 const std::vector<double>& aVal,
 int nStrideVal, int nOffset)
 {
 std::ofstream fout(fpath);
 fout << "# vtk DataFile Version 2.0" << std::endl;
 fout << "hexahedron" << std::endl;
 fout << "ASCII" << std::endl;
 fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
 const int np = aXYZ.size()/3;
 fout << "POINTS " << np << " float" << std::endl;
 for(int ip=0;ip<np;++ip){
 fout << aXYZ[ip*3+0] << " " << aXYZ[ip*3+1] << " " << aXYZ[ip*3+2] << std::endl;
 }
 const int nnoel = 8;
 const int nelem = aHex.size()/nnoel;
 fout << "CELLS " << nelem << " " << nelem*(nnoel+1) << std::endl;
 for(int ielem=0;ielem<nelem;++ielem){
 fout << nnoel << " ";
 fout << aHex[ielem*8+0] << " " << aHex[ielem*8+1] << " " << aHex[ielem*8+3] << " " << aHex[ielem*8+2] << " ";
 fout << aHex[ielem*8+4] << " " << aHex[ielem*8+5] << " " << aHex[ielem*8+7] << " " << aHex[ielem*8+6] << std::endl;
 }
 fout << "CELL_TYPES " << nelem << std::endl;
 for(int ih=0;ih<nelem;++ih){
 fout << "12" << std::endl;
 }
 }
 
 void WriteVTK_TetVec
 (const std::string& fpath,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTet,
 const std::vector<double>& aVal,
 int nStrideVal, int nOffset)
 {
 std::ofstream fout(fpath);
 fout << "# vtk DataFile Version 2.0" << std::endl;
 fout << "velocity" << std::endl;
 fout << "ASCII" << std::endl;
 fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
 const int np = aXYZ.size()/3;
 fout << "POINTS " << np << " float" << std::endl;
 for(int ihp=0;ihp<np;++ihp){
 fout << aXYZ[ihp*3+0] << " " << aXYZ[ihp*3+1] << " " << aXYZ[ihp*3+2] << std::endl;
 }
 const int nnoel = 4;
 const int nelem = aTet.size()/nnoel;
 fout << "CELLS " << nelem << " " << nelem*(nnoel+1) << std::endl;
 for(int ielem=0;ielem<nelem;++ielem){
 fout << nnoel;
 for(int inoel=0;inoel<nnoel;++inoel){
 fout << " " << aTet[ielem*nnoel+inoel];
 }
 fout << std::endl;
 }
 fout << "CELL_TYPES " << nelem << std::endl;
 for(int ielem=0;ielem<nelem;++ielem){
 fout << "10" << std::endl;
 }
 fout << "POINT_DATA " << np << std::endl;
 fout << "VECTORS point_vectors float" << std::endl;
 for(int ip=0;ip<np;++ip){
 fout << aVal[ip*nStrideVal+nOffset+0] << " ";
 fout << aVal[ip*nStrideVal+nOffset+1] << " ";
 fout << aVal[ip*nStrideVal+nOffset+2] << std::endl;
 }
 }
 
 void WriteVTK_TriScalar
 (const std::string& fpath,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<double>& aVal,
 int nStrideVal, int nOffset)
 {
 std::ofstream fout(fpath);
 fout << "# vtk DataFile Version 2.0" << std::endl;
 fout << "hexahedron" << std::endl;
 fout << "ASCII" << std::endl;
 fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
 const int np = aXYZ.size()/3;
 fout << "POINTS " << np << " float" << std::endl;
 for(int ip=0;ip<np;++ip){
 fout << aXYZ[ip*3+0] << " " << aXYZ[ip*3+1] << " " << aXYZ[ip*3+2] << std::endl;
 }
 const int nnoel = 3;
 const int nelem = aTri.size()/nnoel;
 fout << "CELLS " << nelem << " " << nelem*(nnoel+1) << std::endl;
 for(int ielem=0;ielem<nelem;++ielem){
 fout << nnoel << " ";
 fout << aTri[ielem*nnoel+0] << " " << aTri[ielem*nnoel+1] << " " << aTri[ielem*nnoel+2] << std::endl;
 }
 fout << "CELL_TYPES " << nelem << std::endl;
 for(int ih=0;ih<nelem;++ih){
 fout << "5" << std::endl;
 }
 fout << "POINT_DATA " << np << std::endl;
 fout << "SCALARS point_scalars float 1" << std::endl;
 fout << "LOOKUP_TABLE default" << std::endl;
 for(int ip=0;ip<np;++ip){
 fout << aVal[ip*nStrideVal+nOffset] << std::endl;
 }
 }
 
 void WriteVTK_QuadScalar
 (const std::string& fpath,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aQuad,
 const std::vector<double>& aVal,
 int nStrideVal, int nOffset)
 {
 std::ofstream fout(fpath);
 fout << "# vtk DataFile Version 2.0" << std::endl;
 fout << "hexahedron" << std::endl;
 fout << "ASCII" << std::endl;
 fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
 const int np = aXYZ.size()/3;
 fout << "POINTS " << np << " float" << std::endl;
 for(int ip=0;ip<np;++ip){
 fout << aXYZ[ip*3+0] << " " << aXYZ[ip*3+1] << " " << aXYZ[ip*3+2] << std::endl;
 }
 const int nnoel = 4;
 const int nelem = aQuad.size()/nnoel;
 fout << "CELLS " << nelem << " " << nelem*(nnoel+1) << std::endl;
 for(int ielem=0;ielem<nelem;++ielem){
 fout << nnoel << " ";
 fout << aQuad[ielem*nnoel+0] << " ";
 fout << aQuad[ielem*nnoel+1] << " ";
 fout << aQuad[ielem*nnoel+2] << " ";
 fout << aQuad[ielem*nnoel+3] << std::endl;
 }
 fout << "CELL_TYPES " << nelem << std::endl;
 for(int ih=0;ih<nelem;++ih){
 fout << "9" << std::endl;
 }
 fout << "POINT_DATA " << np << std::endl;
 fout << "SCALARS point_scalars float 1" << std::endl;
 fout << "LOOKUP_TABLE default" << std::endl;
 for(int ip=0;ip<np;++ip){
 fout << aVal[ip*nStrideVal+nOffset] << std::endl;
 }
 }
 */

DFM2_INLINE void delfem2::WriteVTK_MapTriScalar(
    const std::string &fpath,
    const std::string &name,
    const std::vector<double> &aXYZ,
    const std::vector<int> &map,
    const std::vector<int> &aTri,
    const std::vector<double> &aVal,
    int nStrideVal,
    int nOffset) {
  std::ofstream fout(fpath.c_str());
  fout << "# vtk DataFile Version 2.0" << std::endl;
  fout << name << std::endl;
  fout << "ASCII" << std::endl;
  fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
  const size_t np = map.size();
  fout << "POINTS " << np << " float" << std::endl;
  for (unsigned int jp = 0; jp < np; ++jp) {
    int ip = map[jp];
    fout << aXYZ[ip * 3 + 0] << " " << aXYZ[ip * 3 + 1] << " " << aXYZ[ip * 3 + 2] << std::endl;
  }
  const int nnoel = 3;
  const std::size_t nelem = aTri.size() / nnoel;
  fout << "CELLS " << nelem << " " << nelem * (nnoel + 1) << std::endl;
  for (unsigned int ielem = 0; ielem < nelem; ++ielem) {
    fout << nnoel << " ";
    fout << aTri[ielem * nnoel + 0] << " " << aTri[ielem * nnoel + 1] << " " << aTri[ielem * nnoel + 2] << std::endl;
  }
  fout << "CELL_TYPES " << nelem << std::endl;
  for (unsigned int ih = 0; ih < nelem; ++ih) {
    fout << "5" << std::endl;
  }
  fout << "POINT_DATA " << np << std::endl;
  fout << "SCALARS point_scalars float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for (unsigned int jp = 0; jp < np; ++jp) {
    int ip = map[jp];
    fout << aVal[ip * nStrideVal + nOffset] << std::endl;
  }
}

DFM2_INLINE void delfem2::ReadVTK(
    std::vector<double> &aXYZ,
    int &ielemtype,
    std::vector<int> &aElem,
    std::vector<double> &aPointVal,
    const std::string &fpath) {
  std::ifstream fin(fpath.c_str());

  const int nbuff = 512;
  char buff[nbuff];
  fin.getline(buff, nbuff);
  fin.getline(buff, nbuff);
  fin.getline(buff, nbuff);
  fin.getline(buff, nbuff);

  std::string str0, str1, str2, str3;
  fin >> str0 >> str1 >> str2;
  //  std::cout << str0 << " " << str1 << " " << str2 << std::endl;

  int np = msh_iomisc::MyStoi_(str1);
  aXYZ.resize(np * 3);
  for (int ip = 0; ip < np; ++ip) {
    fin >> str0 >> str1 >> str2;
    aXYZ[ip * 3 + 0] = msh_iomisc::MyStod_(str0);
    aXYZ[ip * 3 + 1] = msh_iomisc::MyStod_(str1);
    aXYZ[ip * 3 + 2] = msh_iomisc::MyStod_(str2);
  }

  fin >> str0 >> str1 >> str2;
  //  std::cout << str0 << " " << str1 << " " << str2 << std::endl;
  int nelem = msh_iomisc::MyStoi_(str1);
  int nbuffelem = msh_iomisc::MyStoi_(str2);
  std::vector<int> aElembuff(nbuffelem);
  for (int i = 0; i < nbuffelem; ++i) { fin >> aElembuff[i]; }

  fin >> str0 >> str1;
  //  std::cout << str0 << " " << str1 << std::endl;
  int nelem1 = msh_iomisc::MyStoi_(str1);
  std::vector<int> aElemType(nelem1);
  for (int i = 0; i < nelem1; ++i) { fin >> aElemType[i]; }

  fin >> str0 >> str1;
  if (!fin.eof()) {
    //  std::cout << str0 << " " << str1 << std::endl;
    int np1 = msh_iomisc::MyStoi_(str1);
    assert(np1 == np);
    fin >> str0;
    if (str0 == "VECTORS") {
      fin >> str1 >> str2;
      aPointVal.resize(np1 * 3);
      for (int i = 0; i < np1 * 3; ++i) { fin >> aPointVal[i]; }
    }
    if (str0 == "SCALARS") {
      fin >> str1 >> str2 >> str3;
      fin >> str1 >> str2;
      aPointVal.resize(np1);
      for (int i = 0; i < np1; ++i) { fin >> aPointVal[i]; }
    }
  } else {
    aPointVal.resize(0);
  }
  fin.close();
  // ---------
  assert(nelem1 == nelem);
  ielemtype = aElemType[0];
  for (int ie = 0; ie < nelem; ++ie) {
    assert(aElemType[ie] == ielemtype);
  }
  int nnoel = 0;
  if (ielemtype == 5) { nnoel = 3; } // Tri
  if (ielemtype == 9) { nnoel = 4; } // Quad
  if (ielemtype == 10) { nnoel = 4; } // Tet
  if (ielemtype == 12) { nnoel = 8; } // Hex
  assert((int) aElembuff.size() == nelem * (nnoel + 1));
  aElem.resize(nelem * nnoel);
  for (int ie = 0; ie < nelem; ++ie) {
    for (int ine = 0; ine < nnoel; ++ine) {
      aElem[ie * nnoel + ine] = aElembuff[ie * (nnoel + 1) + ine + 1];
    }
  }
}

DFM2_INLINE void delfem2::Read_MeshTri3D_Nas(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    const char *path) {
  aXYZ.clear();
  aTri.clear();
  /////
//  FILE *fp = fopen(path,"r");
  std::ifstream fin(path);
  if (!fin.is_open()) {
    printf("fail to open file");
    return;
  }
  // -----
  std::vector<int> map01;
  const int nbuff = 256;
  std::vector<char> buff(nbuff);
  char buff1[16], buff2[16], buff3[16], buff4[16], buff5[16];
  while (true) {
    fin.getline(buff.data(), nbuff);
    if (fin.eof()) { break; } // eof
    if (buff[0] == 'G' && buff[1] == 'R' && buff[2] == 'I' && buff[3] == 'D') {
      for (int i = 0; i < 8; ++i) { buff1[i] = buff[i + 8]; }
      buff1[8] = '\0';
      for (int i = 0; i < 8; ++i) { buff2[i] = buff[i + 24]; }
      buff2[8] = '\0';
      for (int i = 0; i < 8; ++i) { buff3[i] = buff[i + 32]; }
      buff3[8] = '\0';
      for (int i = 0; i < 8; ++i) { buff4[i] = buff[i + 40]; }
      buff4[8] = '\0';
      const int i0 = atoi(buff1) - 1;
      const double f2 = atof(buff2);
      const double f3 = atof(buff3);
      const double f4 = atof(buff4);
      map01.push_back(i0);
      aXYZ.push_back(f2);
      aXYZ.push_back(f3);
      aXYZ.push_back(f4);
      continue;
    }
    if (buff[0] == 'C' && buff[1] == 'T' && buff[2] == 'R' && buff[3] == 'I' && buff[4] == 'A' && buff[5] == '3') {
      for (int i = 0; i < 8; ++i) { buff1[i] = buff[i + 8]; }
      buff1[8] = '\0';
      for (int i = 0; i < 8; ++i) { buff2[i] = buff[i + 16]; }
      buff2[8] = '\0';
      for (int i = 0; i < 8; ++i) { buff3[i] = buff[i + 24]; }
      buff3[8] = '\0';
      for (int i = 0; i < 8; ++i) { buff4[i] = buff[i + 32]; }
      buff4[8] = '\0';
      for (int i = 0; i < 8; ++i) { buff5[i] = buff[i + 40]; }
      buff5[8] = '\0';
      const int i3 = atoi(buff3) - 1;
      const int i4 = atoi(buff4) - 1;
      const int i5 = atoi(buff5) - 1;
      aTri.push_back(i3);
      aTri.push_back(i4);
      aTri.push_back(i5);
      continue;
    }
  }
  int np1 = map01[0];
  for (int ip1 : map01) {
    np1 = (ip1 > np1) ? ip1 : np1;
  }
  np1 += 1;
  std::vector<unsigned int> map10(np1, UINT_MAX);
  for (unsigned int ip0 = 0; ip0 < map01.size(); ++ip0) {
    unsigned int ip1 = map01[ip0];
    map10[ip1] = ip0;
  }
  for (unsigned int &i : aTri) {
    unsigned int j1 = i;
    i = map10[j1];
  }
}

