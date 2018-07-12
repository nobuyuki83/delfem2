#ifndef meshio_hpp
#define meshio_hpp

#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////
// IO functions

void saveMesh(std::ofstream& fout,
              const std::vector<double>& aXYZ,
              const std::vector<int>& aTri,
              int ndim=3);
void loadMesh(std::ifstream& fin,
              std::vector<double>& aXYZ,
              std::vector<int>& aTri);

/////////////
// TetGen

void Load_TetMesh_TetGen(const std::string& fname,
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
void Load_Ply(const std::string& fname,
              std::vector<double>& aXYZ,
              std::vector<int>& aTri);
void Write_Ply(const std::string& fname,
               unsigned int nXYZ, double* paXYZ,
               unsigned int nTri, unsigned int* paTri);
void Load_Ply(const std::string& fname,
              unsigned int& nnode_, double*& pXYZs_,
              unsigned int& ntri_,  unsigned int*& aTri_);

/////////////
// Obj

void Write_Obj(const std::string& str,
               const std::vector<double>& aXYZ,
               const std::vector<int>& aTri);
void Write_Obj_Quad(const std::string& str,
                    const std::vector<double>& aXYZ,
                    const std::vector<int>& aQuad);
void Write_Obj(const std::string& str, // mixed elem
               const std::vector<double>& aXYZ,
               const std::vector<int>& aElemInd,
               const std::vector<int>& aElem);
void Write_Obj(const std::string& str,
               const std::vector< std::pair< std::vector<double>, std::vector<int> > >& aMesh);
void Load_Obj(const std::string& fname,
              std::vector<double>& aXYZ,
              std::vector<int>& aTri);
void Load_Obj_Quad(const std::string& fname,
                   std::vector<double>& aXYZ,
                   std::vector<int>& aQuad);
void Load_Obj2(const std::string& fname,
               std::vector<double>& aXYZ,
               std::vector<int>& aTri);
void Load_Obj(std::stringstream& ssobj,
              std::vector<double>& aXYZ,
              std::vector<int>& aTri);

class CTriGroup{
public:
  std::string name_group;
  std::string name_mtl;
  int imtl;
  std::vector<int> aTriVtx;
  std::vector<int> aTriNrm;
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
                     const std::vector<double>& aXYZ);

// 5 : VTK_TRIANGLE
// 9 : VTK_QUAD
// 10: VTK_TETRA
// 12: VTK_HEXAHEDRON
// 13: VTK_WEDGE
// 14: VTK_PYRAMD
void WriteVTK_Cells(std::ofstream& fout,
                    int elem_type,
                    const std::vector<int>& aElem);
void WriteVTK_Cells(std::ofstream& fout,
                    const std::vector<int>& aTet,
                    const std::vector<int>& aPyrm,
                    const std::vector<int>& aPrsm);
void WriteVTK_Data_PointVec(std::ofstream& fout,
                            int np,
                            const std::vector<double>& aVal,
                            int nStrideVal=3, int nOffset=0);
void WriteVTK_Data_PointScalar(std::ofstream& fout,
                               int np,
                               const std::vector<double>& aVal,
                               int nStrideVal=1, int nOffset=0);
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

#endif /* meshio_hpp */
