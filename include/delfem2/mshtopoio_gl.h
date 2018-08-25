#ifndef MSHTOPOIO_GL_H
#define MSHTOPOIO_GL_H

#include <string>
#include <vector>
#include <fstream>

#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/funcs_gl.h"

class CMeshElem{
public:
  CMeshElem(){}
  CMeshElem(const std::string& fpath){
    this->Read(fpath);
  }
  void Read(const std::string& fname){
    Read_Ply(fname, aPos, aElem);
    elem_type = MESHELEM_TRI;
    ndim = 3;
  }
  void DrawFace_ElemWiseNorm(){
    if( elem_type == MESHELEM_TRI ){
      if( ndim == 3 ){ DrawMeshTri3D_FaceNorm(aPos, aElem); }
    }
    else if( elem_type == MESHELEM_QUAD ){
      if( ndim == 3 ){ DrawMeshQuad3D_FaceNorm(aPos, aElem); }
    }
  }
  void DrawEdge(){
    if( elem_type == MESHELEM_TRI ){
      if( ndim == 3 ){ DrawMeshTri3D_Edge(aPos, aElem); }
    }
    else if( elem_type == MESHELEM_QUAD ){
      if( ndim == 3 ){ DrawMeshQuad3D_Edge(aPos, aElem); }
    }
  }
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
};

#endif
