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
#include "delfem2/vec3.h"
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
  void InitGL();
  std::vector<double> AABBVec3() const {
    std::vector<double> mm(6);
    mm[0] = +1;
    mm[1] = -1;
    return mm;
  }
  void SetZeroToDepth(){ for(unsigned int i=0;i<aZ.size();++i){ aZ[i] = 0.0; } }
  void GetMVPG(double mMVPG[16]) const {
    double mMVP[16]; MatMat4(mMVP, mMV,mP);
    const double tmp0 = nResX*0.5;
    const double tmp1 = nResY*0.5;
    double mG[16] = {
        tmp0,    0,   0, 0,
        0,    tmp1,   0, 0,
        0,       0, 0.5, 0,
        tmp0-0.5, tmp1-0.5, 0.5, 1 };
    MatMat4(mMVPG, mMVP,mG);
  }
  /**
  * @brief update the bounding box by adding points
  * @param pmin (in/out) lower coner
  * @param pmax (in/out) upper corner
  * @details if( pmin[0] > pmax[0] ) this bounding box is empty
  */
  void BoundingBox3(double* pmin, double* pmax) const;
  // ----------------------
  void SaveDepthCSV(
      const std::string& path) const;

  void SetTextureProperty(unsigned int nw, unsigned int nh, bool is_rgba_8ui_)
  {
    this->nResX = nw;
    this->nResY = nh;
    this->is_rgba_8ui = is_rgba_8ui_;
  }
  void Start();
  void End();
  void ExtractFromTexture_Depth(std::vector<float>& aZ);
  void ExtractFromTexture_RGBA8UI(std::vector<std::uint8_t>& aRGBA);
  void ExtractFromTexture_RGBA32F(std::vector<float>& aRGBA);
  void GetDepth()
  {
    CRender2Tex::ExtractFromTexture_Depth(aZ);
  }

  void GetColor()
  {
    if( is_rgba_8ui ){
      CRender2Tex::ExtractFromTexture_RGBA8UI(aRGBA_8ui);
    }
    else{
      CRender2Tex::ExtractFromTexture_RGBA32F(aRGBA_32f);
    }
  }
public:
  unsigned int nResX;
  unsigned int nResY;
  // ------------------
  unsigned int id_tex_color;
  unsigned int id_tex_depth;
  unsigned int id_framebuffer;
  // --------
  std::vector<float> aZ;
  bool is_rgba_8ui;
  std::vector<unsigned char> aRGBA_8ui;
  std::vector<float> aRGBA_32f;
  //
  double mMV[16]; // affine matrix
  double mP[16]; // affine matrix
protected:
  int view[4]; // viewport information
};

/**
 * @brief project input point to the depth surface
 * @param[in] ps the point to project
 */
bool GetProjectedPoint(
                       CVec3d& p0,
                       CVec3d& n0,
                       const CVec3d& ps,
                       const CRender2Tex& smplr);


} // opengl
} // delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/r2t_gl.cpp"
#endif

#endif /* depth_hpp */
