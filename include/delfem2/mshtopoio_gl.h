#ifndef MSHTOPOIO_GL_H
#define MSHTOPOIO_GL_H

#include <string>
#include <vector>
#include <fstream>

#if defined(__APPLE__) && defined(__MACH__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else // linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/bv.h"

class CMeshElem{
public:
  CMeshElem(){}
  CMeshElem(const std::string& fpath){
    this->Read(fpath);
  }
  void Draw(){
    glEnable(GL_LIGHTING);
    this->DrawFace_ElemWiseNorm();
    glDisable(GL_LIGHTING);
    glLineWidth(1);
    this->DrawEdge();
  }
  CBV3D_AABB AABB3() const{
    CBV3D_AABB aabb;
    double cw[6]; GetCenterWidth(cw, aPos);
    aabb.x_min = cw[0]-0.5*cw[3];
    aabb.x_max = cw[0]+0.5*cw[3];
    aabb.y_min = cw[1]-0.5*cw[4];
    aabb.y_max = cw[1]+0.5*cw[4];
    aabb.z_min = cw[2]-0.5*cw[5];
    aabb.z_max = cw[2]+0.5*cw[5];
    return aabb;
  }
  void Read(const std::string& fname){
    Read_Ply(fname, aPos, aElem);
    elem_type = MESHELEM_TRI;
    ndim = 3;
  }
  void DrawFace_ElemWiseNorm(){
    if( elem_type == MESHELEM_TRI ){
      if( ndim == 3 ){ DrawMeshTri3D_FaceNorm(aPos, aElem); }
      if( ndim == 2 ){ DrawMeshTri2D_Face(aElem,aPos); }
    }
    else if( elem_type == MESHELEM_QUAD ){
      if( ndim == 3 ){ DrawMeshQuad3D_FaceNorm(aPos, aElem); }
    }
  }
  void DrawEdge(){
    if( elem_type == MESHELEM_TRI ){
      if( ndim == 3 ){ DrawMeshTri3D_Edge(aPos, aElem); }
      if( ndim == 2 ){ DrawMeshTri2D_Edge(aElem, aPos); }
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
