/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// TODO: having MVP matrix instead of axes

#ifndef DFM2_RENDER2TEX_GL_H
#define DFM2_RENDER2TEX_GL_H

#include "delfem2/mat4.h"
#include "delfem2/dfm2_inline.h"
#include <stdio.h>
#include <vector>

namespace delfem2 {

void Mat4_OrthongoalProjection_AffineTrans(
    double mMV[16],
    double mP[16],
    const double origin[3],
    const double az[3],
    const double ax[3],
    unsigned int nResX,
    unsigned int nResY,
    double lengrid,
    double z_range);

namespace opengl {

class CRender2Tex
{
public:
  CRender2Tex() :
      mMV{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1},
      mP{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}
  {
    nResX=0;
    nResY=0;
    is_rgba_8ui=false;
    id_tex_color = 0;
    id_tex_depth = 0;
  }
  // --------------------------
  virtual void InitGL();
  std::vector<double> AABBVec3() const {
    std::vector<double> mm(6);
    mm[0] = +1;
    mm[1] = -1;
    return mm;
  }
  void GetMVPG(double mMVPG[16]) const {
    double mMVP[16]; MatMat4(mMVP, mMV,mP);
    double tmp0 = nResX*0.5;
    double tmp1 = nResY*0.5;
    double mG[16] = {
        tmp0,    0,   0, 0,
        0,    tmp1,   0, 0,
        0,       0, 0.5, 0,
        tmp0-0.5, tmp1-0.5, 0.5, 1 };
    MatMat4(mMVPG, mMVP,mG);
//    p0.p[0] = (p0.p[0]+1)*sampler.nResX*0.5 - 0.5;
//    p0.p[1] = (p0.p[1]+1)*sampler.nResY*0.5 - 0.5;
//    p0.p[2] = (p0.p[2]+1)*0.5;
  }
  // ----------------------
  /*
  void AffMatT3f_MVP(
      float mMV[16],
      float p[16]) const;
      */
  void SaveDepthCSV(
      const std::string& path) const;
  /*
  void SetCoord(
      double elen, double depth_max,
      const std::vector<double>& orgPrj,
      const std::vector<double>& dirPrj,
      const std::vector<double>& dirWidth);
      */
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
  //
  double mMV[16]; // affine matrix
  double mP[16]; // affine matrix
protected:
  int view[4]; // viewport information
};


} // opengl
} // delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/r2t_gl.cpp"
#endif

#endif /* depth_hpp */
