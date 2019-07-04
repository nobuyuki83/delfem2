/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef MSHTOPOIO_H
#define MSHTOPOIO_H

#include <string>
#include <vector>
#include <fstream>

#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/funcs.h"

void MeshTri3D_GeodesicPolyhedron(std::vector<double>& aXYZ1,
                                  std::vector<unsigned int>& aTri1);

class CMeshElem{
public:
  CMeshElem(){
    color_face.resize(4);
    color_face[0] = 0.8;
    color_face[1] = 0.8;
    color_face[2] = 0.8;
    color_face[3] = 1.0;
    is_draw_edge = true;
  }
  CMeshElem(const std::string& fpath){
    color_face.resize(4);
    color_face[0] = 0.8;
    color_face[1] = 0.8;
    color_face[2] = 0.8;
    color_face[3] = 1.0;
    is_draw_edge = true;
    this->Read(fpath);
  }
//  void Draw() const;
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
    std::string sExt = pathGetExtension(fname);
    if( sExt == "ply") {
      Read_Ply(fname, aPos, aElem);
    }
    else if( sExt == "obj") {
      Read_Obj(fname, aPos, aElem);
    }
    elem_type = MESHELEM_TRI;
    ndim = 3;
  }
  void Write_Obj(const std::string& fname){
    ::Write_Obj(fname,
                aPos.data(), aPos.size()/3,
                aElem.data(), aElem.size()/3);
  }
//  void DrawFace_ElemWiseNorm() const;
//  void DrawEdge() const;
  CMeshElem Subdiv(){
    CMeshElem em;
    if( elem_type == MESHELEM_QUAD ){
      const std::vector<double>& aXYZ0 = this->aPos;
      const std::vector<unsigned int>& aQuad0 = this->aElem;
      em.elem_type = MESHELEM_QUAD;
      em.ndim = 3;
      std::vector<unsigned int>& aQuad1 = em.aElem;
      std::vector<int> aEdgeFace0;
      std::vector<int> psupIndQuad0, psupQuad0;
      QuadSubdiv(aQuad1,
                 psupIndQuad0,psupQuad0, aEdgeFace0,
                 aQuad0.data(), aQuad0.size()/4, aXYZ0.size()/3);
      ///////
      std::vector<double>& aXYZ1 = em.aPos;
      SubdivisionPoints_QuadCatmullClark(aXYZ1,
                                         aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,
                                         aQuad0.data(),aQuad0.size()/4,
                                         aXYZ0.data(),aXYZ0.size()/3);
    }
    return em;
  }
  void ScaleXYZ(double s){
    Scale(s,aPos);
  }
public:
  MESHELEM_TYPE elem_type;
  std::vector<unsigned int> aElem;
  /////
  int ndim;
  std::vector<double> aPos;
  /////
  std::vector<float> color_face;
  bool is_draw_edge;
};

class CMaterial{
public:
  std::string name_mtl;
  float Kd[4];
  float Ka[4];
  float Ks[4];
  float Ke[4];
  float Ns;
  int illum;
  std::string map_Kd;
public:
  void GL() const;
};

class CMeshMultiElem{
public:
  void ReadObj(const std::string& fname);
//  void Draw() const;
  std::vector<double> AABB3_MinMax() const;
  void ScaleXYZ(double s);
  void TranslateXYZ(double x, double y, double z);
public:
  std::vector<double> aXYZ;
  std::vector<double> aNorm;
  std::vector<CTriGroup> aObjGroupTri;
  std::vector<CMaterial> aMaterial;
};


void Load_Mtl(const std::string& fname,
              std::vector<CMaterial>& aMtl);

#endif
