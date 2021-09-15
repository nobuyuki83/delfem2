/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_MSHTOPOIO_H
#define DFM2_MSHTOPOIO_H

#include <string>
#include <vector>
#include <fstream>

#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_io_misc.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshsubdiv.h"
#include "delfem2/jagarray.h"
#include "delfem2/points.h"
#include "delfem2/str.h"
#include "delfem2/file.h"

void MeshTri3D_GeodesicPolyhedron(
    std::vector<double>& aXYZ1,
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
    double c[3], w[3];
    delfem2::CenterWidth_Points3(c, w,
                                 aPos);
    std::vector<double> aabb(6);
    aabb[0] = c[0]-0.5*w[0];
    aabb[1] = c[0]+0.5*w[0];
    aabb[2] = c[1]-0.5*w[1];
    aabb[3] = c[1]+0.5*w[1];
    aabb[4] = c[2]-0.5*w[2];
    aabb[5] = c[2]+0.5*w[2];
    return aabb;
  }
  void Read(const std::string& fname){
    std::string sExt = delfem2::pathGetExtension(fname);
    if( sExt == "ply") {
      delfem2::Read_Ply(
          aPos, aElem,
          std::filesystem::path(fname));
    }
    else if( sExt == "obj") {
      delfem2::Read_Obj(fname, aPos, aElem);
    }
    elem_type = delfem2::MESHELEM_TRI;
    ndim = 3;
  }
  void Write_Obj(const std::string& fname){
    delfem2::Write_Obj(fname,
                aPos.data(), aPos.size()/3,
                aElem.data(), aElem.size()/3);
  }
//  void DrawFace_ElemWiseNorm() const;
//  void DrawEdge() const;
  CMeshElem Subdiv(){
    CMeshElem em;
    if( elem_type == delfem2::MESHELEM_QUAD ){
      const std::vector<double>& aXYZ0 = this->aPos;
      const std::vector<unsigned int>& aQuad0 = this->aElem;
      em.elem_type = delfem2::MESHELEM_QUAD;
      em.ndim = 3;
      std::vector<unsigned int>& aQuad1 = em.aElem;
      std::vector<unsigned int> aEdgeFace0;
      std::vector<unsigned int> psupIndQuad0, psupQuad0;
      delfem2::SubdivTopo_MeshQuad(
          aQuad1,
          psupIndQuad0,psupQuad0, aEdgeFace0,
          aQuad0.data(), aQuad0.size()/4, aXYZ0.size()/3);
      // ----------------
      std::vector<double>& aXYZ1 = em.aPos;
      delfem2::SubdivisionPoints_QuadCatmullClark(
          aXYZ1,
          aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,
          aQuad0.data(),aQuad0.size()/4,
          aXYZ0.data(),aXYZ0.size()/3);
    }
    return em;
  }
  void ScaleXYZ(double s){
    delfem2::Scale_PointsX(aPos,
                           s);
  }
public:
  delfem2::MESHELEM_TYPE elem_type;
  std::vector<unsigned int> aElem;
  //
  int ndim;
  std::vector<double> aPos;
  //
  std::vector<float> color_face;
  bool is_draw_edge;
};



#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mshtopoio.cpp"
#endif

#endif
