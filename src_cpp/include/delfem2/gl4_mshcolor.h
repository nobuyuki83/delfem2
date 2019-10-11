/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef GL4_MSH_H
#define GL4_MSH_H

#include <stdio.h>
#include <vector>

#include "delfem2/color.h"
#include "delfem2/gl4_funcs.h" // CGL4_VAO_Mesh

class CShader_TriMesh{
public:
  void Initialize(std::vector<double>& aXYZd,
                  std::vector<unsigned int>& aTri);
  void UpdateVertex(std::vector<double>& aXYZd,
                    std::vector<unsigned int>& aTri);
  void Compile();
  void Draw(float mP[16], float mMV[16]);
  
public:
  CGL4_VAO_Mesh vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color;
};


class CShader_TriMesh_Scalar{
public:
  CShader_TriMesh_Scalar(){
    val_min = 0.0;
    val_max = 1.0;
    color_min = CColor::Gray(0.0);
    color_max = CColor::Gray(1.0);
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
  CColor color_min, color_max;
};


#endif /* gl4_msh_hpp */
