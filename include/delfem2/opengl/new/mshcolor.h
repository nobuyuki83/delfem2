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
  void Initialize(std::vector<double>& aXYZd);
  void UpdateVertex(std::vector<double>& aXYZd);
  void Draw(float mP[16], float mMV[16]) const;
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

  void Initialize(std::vector<double>& aXYZd,
                  std::vector<unsigned int>& aLine);

  void UpdateVertex(std::vector<double>& aXYZd,
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

  void Initialize(std::vector<double>& aXYZd,
                  std::vector<unsigned int>& aTri);
  void UpdateVertex(std::vector<double>& aXYZd,
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
  void Initialize(std::vector<double>& aPosD,
                  unsigned int ndim,
                  std::vector<unsigned int>& aTri,
                  std::vector<double>& aValD);
  void UpdateVertex(std::vector<double>& aPosD,
                    unsigned int ndim,
                    std::vector<double>& aValD);
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
  void Initialize(std::vector<double>& aPosD,
                  unsigned int ndim,
                  std::vector<unsigned int>& aTri,
                  std::vector<double>& aDispD);
  void UpdateVertex(std::vector<double>& aPosD,
                    unsigned int ndim,
                    std::vector<double>& aDispD);
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
