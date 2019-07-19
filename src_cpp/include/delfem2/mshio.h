/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef MSHIO_H
#define MSHIO_H

#include <stdio.h>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////
// IO functions

void Write_MeshTri3D(std::ofstream& fout,
                    const std::vector<double>& aXYZ,
                    const std::vector<int>& aTri,
              int ndim=3);
void Read_MeshTri3D(std::ifstream& fin,
                    std::vector<double>& aXYZ,
                    std::vector<int>& aTri);

/////////////
// TetGen

void Read_MeshTet3D_TetGen(const std::string& fname,
                           std::vector<double>& aXYZ,
                           std::vector<int>& aTet,
                           std::vector<int>& aTri);

/////////////
// STL

void Write_STL(const std::string& str,
               const std::vector<double>& aXYZ,
               const std::vector<int>& aTri);


/////////////
// PLY

void Write_Ply(const std::string& fname,
               const std::vector<double>& aXYZ,
               const std::vector<int>& aTri);
void Write_Ply(const std::string& fname,
               unsigned int nXYZ, double* paXYZ,
               unsigned int nTri, unsigned int* paTri);
void Read_Ply(const std::string& fname,
              int& nnode_, double*& pXYZs_,
              int& ntri_,  unsigned int*& aTri_);
void Read_Ply(const std::string& fname,
              std::vector<double>& aXYZ,
              std::vector<unsigned int>& aTri);

/////////////
// Obj

void Write_Obj(const std::string& str,
               const double* aXYZ, int nXYZ,
               const unsigned int* aTri, int nTri);
void Write_Obj_Quad(const std::string& str,
                    const std::vector<double>& aXYZ,
                    const std::vector<int>& aQuad);
void Write_Obj(const std::string& str, // mixed elem
               const std::vector<double>& aXYZ,
               const std::vector<int>& aElemInd,
               const std::vector<int>& aElem);
void Write_Obj(const std::string& str,
               const std::vector< std::pair< std::vector<double>, std::vector<int> > >& aMesh);
void Read_Obj(const std::string& fname,
              std::vector<double>& aXYZ,
              std::vector<unsigned int>& aTri);
void Read_Obj_Quad(const std::string& fname,
                   std::vector<double>& aXYZ,
                   std::vector<int>& aQuad);
void Read_Obj(std::stringstream& ssobj,
              std::vector<double>& aXYZ,
              std::vector<int>& aTri);
void Read_Obj2(const std::string& fname,
               std::vector<double>& aXYZ,
               std::vector<int>& aTri);
void Read_Obj3(const std::string& fname,
               std::vector<double>& aXYZ,
               std::vector<unsigned int>& aTri);

class CTriGroup{
public:
  std::string name_group;
  std::string name_mtl;
  int imtl;
  std::vector<unsigned int> aTriVtx;
  std::vector<unsigned int> aTriNrm;
};

void Load_Obj(const std::string& fname,
              std::string& fname_mtl,
              std::vector<double>& aXYZ,
              std::vector<double>& aNorm,
              std::vector<CTriGroup>& aTriGroup);

/////////////
// VTK

void WriteVTK_Points(std::ofstream& fout,
                     const std::string& name,
                     const double* pXYZ,
                     int nXYZ,
                     int nDim);



// 5 : VTK_TRIANGLE
// 9 : VTK_QUAD
// 10: VTK_TETRA
// 12: VTK_HEXAHEDRON
// 13: VTK_WEDGE
// 14: VTK_PYRAMD
void WriteVTK_Cells(std::ofstream& fout,
                    int vtk_elem_type,
                    const int* aElem,
                    const int nElem);
void WriteVTK_Cells(std::ofstream& fout,
                    const std::vector<int>& aTet,
                    const std::vector<int>& aPyrm,
                    const std::vector<int>& aPrsm);
void WriteVTK_Data_PointVec(std::ofstream& fout,
                            const double* aVal,
                            int np,
                            int nStrideVal,
                            int ndim);
void WriteVTK_Data_PointScalar(std::ofstream& fout,
                               const double* aVal,
                               int np,
                               int nStrideVal=1);
void WriteVTK_MapTriScalar(const std::string& fpath,
                           const std::string& name,
                           const std::vector<double>& aXYZ,
                           const std::vector<int>& map,
                           const std::vector<int>& aTri,
                           const std::vector<double>& aVal,
                           int nStrideVal, int nOffset);
void ReadVTK(std::vector<double>& aXYZ,
             int& type_elem,
             std::vector<int>& aElem,
             std::vector<double>& aPointVal,
             const std::string& fpath);

void Read_MeshTri3D_Nas(std::vector<double>& aXYZ,
                        std::vector<unsigned int>& aTri,
                        const char* path);

#endif /* meshio_hpp */
