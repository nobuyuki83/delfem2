/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_RENDER2TEX_GL_H
#define DFM2_RENDER2TEX_GL_H

#include "delfem2/dfm2_inline.h"
#include <stdio.h>
#include <vector>

namespace delfem2 {
namespace opengl {

class CRender2Tex
{
public:
  CRender2Tex()
  {
    nResX=0;
    nResY=0;
    is_rgba_8ui=false;
    id_tex_color = 0;
    id_tex_depth = 0;
    lengrid = 0.01;
    origin[0]=0; origin[1]=0; origin[2]=0;
    z_axis[0]=0; z_axis[1]=0; z_axis[2]=1;
    x_axis[0]=1; x_axis[1]=0; x_axis[2]=0;
  }
  // --------------------------
  virtual void InitGL();
  std::vector<double> AABBVec3() const {
    std::vector<double> mm(6);
    mm[0] = +1;
    mm[1] = -1;
    return mm;
  }
  // ----------------------
  void AffMatT3f_MVP(
      float mMV[16],
      float p[16]) const;
  void SaveDepthCSV(
      const std::string& path) const;
  void SetCoord(
      double elen, double depth_max,
      const std::vector<double>& orgPrj,
      const std::vector<double>& dirPrj,
      const std::vector<double>& dirWidth);
  void SetTextureProperty(unsigned int nw, unsigned int nh, bool is_rgba_8ui_)
  {
    this->nResX = nw;
    this->nResY = nh;
    this->is_rgba_8ui = is_rgba_8ui_;
  }
  virtual void Start();
  void End();
  void ExtractFromTexture_Depth(std::vector<float>& aZ);
  void ExtractFromTexture_RGBA8UI(std::vector<std::uint8_t>& aRGBA);
  void ExtractFromTexture_RGBA32F(std::vector<float>& aRGBA);
public:
  unsigned int nResX;
  unsigned int nResY;
  bool is_rgba_8ui;
  // ------------------
  unsigned int id_tex_color;
  unsigned int id_tex_depth;
  unsigned int id_framebuffer;
  double lengrid;
  double z_range;
  double z_axis[3];
  double x_axis[3];
  double origin[3];
  // --------------
protected:
  int view[4]; // viewport information
};


} // opengl
} // delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/r2t_gl.cpp"
#endif

#endif /* depth_hpp */
