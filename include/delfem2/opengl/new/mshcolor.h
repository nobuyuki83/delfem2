/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_NEW_MSHCOLOR_H
#define DFM2_OPENGL_NEW_MSHCOLOR_H

#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/color.h"
#include "delfem2/dfm2_inline.h"
#include <stdio.h>
#include <vector>

// -------------------------------------

namespace delfem2 {
namespace opengl {

class CShader_Points{
public:
  void Compile();

  template <typename REAL>
  void Initialize(
      std::vector<REAL>& aXYZd);

  template <typename REAL>
  void UpdateVertex(
      std::vector<REAL>& aXYZd);

  void Draw(
      float mP[16],
      float mMV[16]) const;
public:
  CGL4_VAO_Mesh vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color;
  unsigned int nPoint = 0;
  delfem2::CColor color_face = delfem2::CColor(0.0,0.0,0.0,0.0);
};

class CShader_LineMesh{
public:
  void Compile();

  template <typename REAL>
  void Initialize(
      std::vector<REAL>& aXYZd,
      std::vector<unsigned int>& aLine);

  template <typename REAL>
  void UpdateVertex(
      std::vector<REAL>& aXYZd,
      std::vector<unsigned int>& aLine);

  void Draw(float mP[16], float mMV[16]) const;

public:
  CGL4_VAO_Mesh vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color;
};

class CShader_TriMesh{
public:
  void Compile();

  template <typename REAL>
  void Initialize(
      std::vector<REAL>& aXYZd,
      std::vector<unsigned int>& aTri);

  template <typename REAL>
  void UpdateVertex(
      std::vector<REAL>& aXYZd,
      std::vector<unsigned int>& aTri);

  void Draw(float mP[16], float mMV[16]) const;
  
public:
  CGL4_VAO_Mesh vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color;
  delfem2::CColor color_face = delfem2::CColor(1.0,0.0,0.0,0.0);
  float line_width = 1.0;
};


class CShader_TriMesh_Scalar{
public:
  CShader_TriMesh_Scalar(){
    val_min = 0.0;
    val_max = 1.0;
    color_min = delfem2::CColor::Gray(0.0);
    color_max = delfem2::CColor::Gray(1.0);
  }

  template <typename REAL>
  void Initialize(
      std::vector<REAL>& aPosD,
      unsigned int ndim,
      std::vector<unsigned int>& aTri,
      std::vector<REAL>& aValD);

  template <typename REAL>
  void UpdateVertex(
      std::vector<REAL>& aPosD,
      unsigned int ndim,
      std::vector<REAL>& aValD);

  void Compile();
  void Draw(float mP[16], float mMV[16]);
  
public:
  CGL4_VAO_Mesh vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color0, Loc_Color1;
  int Loc_ValMin, Loc_ValMax;
  
  double val_min, val_max;
  delfem2::CColor color_min, color_max;
};

class CShader_TriMesh_Disp{
public:
  CShader_TriMesh_Disp(){
  }

  template <typename REAL>
  void Initialize(
      std::vector<REAL>& aPosD,
      unsigned int ndim,
      std::vector<unsigned int>& aTri,
      std::vector<REAL>& aDispD);

  template <typename REAL>
  void UpdateVertex(
      std::vector<REAL>& aPosD,
      unsigned int ndim,
      std::vector<REAL>& aDispD);

  void Compile();
  void Draw(float mP[16], float mMV[16]);
  
public:
  CGL4_VAO_Mesh vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color0, Loc_Color1;
  int Loc_ValMin, Loc_ValMax;
};



}
}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/new/mshcolor.cpp"
#endif

#endif /* gl4_msh_hpp */
