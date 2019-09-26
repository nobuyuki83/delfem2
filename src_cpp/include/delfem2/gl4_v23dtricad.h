/**
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief classes and functions that depend on cad,dtri,vec2,vec3 and OpenGL ver. 4
 */


#ifndef gl4_v23dtricad_h
#define gl4_v23dtricad_h

#include "gl4_funcs.h" // for CGL4_VAO_Mesh

class CCad2D;

class CShader_CCad2D
{
public:
  void MakeBuffer(const CCad2D& cad);
  void Compile(){
    this->Compile_Face();
    this->Compile_Edge();
  }
  void Draw(float mP[16], float mMV[16]) const;
private:
  void Compile_Face();
  void Compile_Edge();
public:
  int shdr0_program; // for face
  int shdr0_Loc_MatrixProjection;
  int shdr0_Loc_MatrixModelView;
  int shdr0_Loc_Color;
  ///
  int shdr1_program; // for face
  int shdr1_Loc_MatrixProjection;
  int shdr1_Loc_MatrixModelView;
  int shdr1_Loc_Color;
  int shdr1_Loc_LineWidth;
  ///
  CGL4_VAO_Mesh vao_face;
  CGL4_VAO_Mesh vao_edge;
};

#endif /* gl4_v23dtricad_h */
