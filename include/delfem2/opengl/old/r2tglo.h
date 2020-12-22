/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file definition of render to texture class CRender2Tex_DrawOldGL
 */

#ifndef DFM2_OPENGL_R2TGLO_GLOLD_H
#define DFM2_OPENGL_R2TGLO_GLOLD_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/opengl/r2t.h"
#include "delfem2/vec3.h"
#include <stdio.h>
#include <vector>
#include <cmath>
#include <cassert>

namespace delfem2 {
namespace opengl {

DFM2_INLINE void SetView(const CRender2Tex& r2t);

class CDrawerOldGL_Render2Tex
{
public:
  CDrawerOldGL_Render2Tex(){}
  ~CDrawerOldGL_Render2Tex(){}
  // ------------
  void Draw(const CRender2Tex& r2t) const;
  void Draw_Axis(const CRender2Tex& r2t) const;
  void Draw_Point(const CRender2Tex& r2t) const;
  void Draw_BoundingBox(const CRender2Tex& r2t) const;
  /**
   * @details before calling this function, bound texture by "glBindTexture" by yourself.
   */
  void Draw_Texture(const CRender2Tex& r2t) const;
  // ------------
  void SetPointColor(double r, double g, double b);
public:
  // -------------------
  double draw_len_axis = 1.0;
  unsigned int pointSize = 3;
  bool isDrawTex = true;
  bool isDrawDepth = true;
  bool isDrawOnlyHitPoints = false;
  std::vector<double> colorPoint = {1,1,1,1};
};



class CRender2Tex_DrawOldGL_BOX
{
public:
  void Draw() const {
    assert( aDrawSampler.size() == aSampler.size() );
    for(unsigned int is=0;is<aSampler.size();++is){
      aDrawSampler[is].Draw(aSampler[is]);
    }
  }
  void Initialize(unsigned int nresX,
                  unsigned int nresY,
                  unsigned int nresZ,
                  double elen);
  
  unsigned int nDivX() const {
    const unsigned int n0 = aSampler[2].nResX;
    assert( aSampler[3].nResX == n0 );
    assert( aSampler[4].nResX == n0 );
    assert( aSampler[5].nResX == n0 );
    return n0;
  }
  unsigned int nDivY() const {
    const unsigned int n0 = aSampler[0].nResX;
    assert( aSampler[1].nResX == n0 );
    assert( aSampler[4].nResY == n0 );
    assert( aSampler[5].nResY == n0 );
    return n0;
  }
  unsigned int nDivZ() const {
    const unsigned int n0 = aSampler[0].nResY;
    assert( aSampler[1].nResY == n0 );
    assert( aSampler[2].nResY == n0 );
    assert( aSampler[3].nResY == n0 );
    return n0;
  }
  double edgeLen() const {
    return this->lengrid;
  }
  void BoundingBox3(double* pmin, double* pmax) const{
    for(const auto& smplr : aSampler ){
      smplr.BoundingBox3(pmin, pmax);
    }
  }
  
public:
  double lengrid;
  std::vector<CDrawerOldGL_Render2Tex> aDrawSampler;
  std::vector<CRender2Tex> aSampler;
};

void CarveVoxelByDepth(
    std::vector<int>& aVal,
    const CRender2Tex_DrawOldGL_BOX& sampler_box);


}
}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/old/r2tglo.cpp"
#endif

#endif /* depth_hpp */
