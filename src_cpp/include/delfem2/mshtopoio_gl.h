#ifndef MSHTOPOIO_GL_H
#define MSHTOPOIO_GL_H

#include <string>
#include <vector>
#include <fstream>

#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
//#include "delfem2/funcs_gl.h"
#include "delfem2/color_gl.h"

void MeshTri3D_GeodesicPolyhedron(std::vector<double>& aXYZ1,
                                  std::vector<int>& aTri1);

class CMeshElem{
public:
  CMeshElem(){
    color_face.resize(4);
    color_face[0] = 0.8;
    color_face[1] = 0.8;
    color_face[2] = 0.8;
    color_face[3] = 1.0;
  }
  CMeshElem(const std::string& fpath){
    color_face.resize(4);
    color_face[0] = 0.8;
    color_face[1] = 0.8;
    color_face[2] = 0.8;
    color_face[3] = 1.0;
    this->Read(fpath);
  }
  void Draw() const;
  std::vector<double> AABB3_MinMax() const{
    double cw[6]; GetCenterWidth(cw, aPos);
    std::vector<double> aabb(6);
    aabb[0] = cw[0]-0.5*cw[3];
    aabb[1] = cw[0]+0.5*cw[3];
    aabb[2] = cw[1]-0.5*cw[4];
    aabb[3] = cw[1]+0.5*cw[4];
    aabb[4] = cw[2]-0.5*cw[5];
    aabb[5] = cw[2]+0.5*cw[5];
    return aabb;
  }
  void Read(const std::string& fname){
    Read_Ply(fname, aPos, aElem);
    elem_type = MESHELEM_TRI;
    ndim = 3;
  }
  void DrawFace_ElemWiseNorm() const;
  void DrawEdge() const;
  CMeshElem Subdiv(){
    CMeshElem em;
    if( elem_type == MESHELEM_QUAD ){
      const std::vector<double>& aXYZ0 = this->aPos;
      const std::vector<int>& aQuad0 = this->aElem;
      em.elem_type = MESHELEM_QUAD;
      em.ndim = 3;
      std::vector<int>& aQuad1 = em.aElem;
      std::vector<int> aEdgeFace0;
      std::vector<int> psupIndQuad0, psupQuad0;
      QuadSubdiv(aQuad1,
                 psupIndQuad0,psupQuad0, aEdgeFace0,
                 aQuad0, aXYZ0.size()/3);
      ///////
      std::vector<double>& aXYZ1 = em.aPos;
      SubdivisionPoints_QuadCatmullClark(aXYZ1,
                                         aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,aQuad0,aXYZ0);
    }
    return em;
  }
  void ScaleXYZ(double s){
    Scale(s,aPos);
  }
public:
  MESHELEM_TYPE elem_type;
  std::vector<int> aElem;
  /////
  int ndim;
  std::vector<double> aPos;
  /////
  std::vector<float> color_face;
};

class CMeshMultiElem{
public:
  void ReadObj(const std::string& fname);
  void Draw() const;
  std::vector<double> AABB3_MinMax() const;
public:
  std::vector<double> aXYZ;
  std::vector<double> aNorm;
  std::vector<CTriGroup> aObjGroupTri;
  std::vector<CMaterial> aMaterial;
};

///////////////////////////////////

class CTriangulationOutput
{
public:
  CMeshElem me;
  std::vector<int> aPtrVtxInd;
  std::vector<int> aVtxInd;
};

CTriangulationOutput Triangulation(const std::vector<double>& aXY,
                                   double edge_length);

#endif
