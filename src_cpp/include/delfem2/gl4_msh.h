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

#include "gl4_funcs.h" // CGL4_VAO_Mesh

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


#endif /* gl4_msh_hpp */
