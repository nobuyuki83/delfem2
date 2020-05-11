/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_RENDER2TEX_GLOLD_H
#define DFM2_RENDER2TEX_GLOLD_H
#include "delfem2/dfm2_inline.h"
#include <stdio.h>
#include <vector>
#include <cmath>
#include "delfem2/opengl/render2tex_gl.h"

namespace delfem2 {
namespace opengl {


class CRender2Tex_DrawOldGL : public CRender2Tex
{
public:
  CRender2Tex_DrawOldGL(){}
  virtual ~CRender2Tex_DrawOldGL(){}
  // ------------
  virtual void InitGL(); // override function
  virtual void Start(); // override function
  // ----------
  void Draw() const;
  void Draw_Axis() const;
  void Draw_Point() const;
  void Draw_BoundingBox() const;
  void SetView();
  // ------------
  void SetPointColor(double r, double g, double b);
  void SetZeroToDepth(){ for(unsigned int i=0;i<aZ.size();++i){ aZ[i] = 0.0; } }
  std::vector<double> getGPos(int ix, int iy) const;
  void GetDepth();
  void GetColor();
public:
  std::vector<float> aZ;
  std::vector<unsigned char> aRGBA_8ui;
  std::vector<float> aRGBA_32f;
  // -------------------
  double draw_len_axis = 1.0;
  unsigned int pointSize = 3;
  bool isDrawTex = true;
  bool isDrawOnlyHitPoints = false;
  std::vector<double> colorPoint = {1,1,1,1};
};
  
class CRender2Tex_DrawOldGL_BOX {
 
public:
  void Draw() const {
    for(auto& smplr: aSampler){
      smplr.Draw();
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
    double l0 = aSampler[0].lengrid;
    assert( fabs(aSampler[1].lengrid-l0) < 1.0e-10 );
    assert( fabs(aSampler[2].lengrid-l0) < 1.0e-10 );
    assert( fabs(aSampler[3].lengrid-l0) < 1.0e-10 );
    assert( fabs(aSampler[4].lengrid-l0) < 1.0e-10 );
    assert( fabs(aSampler[5].lengrid-l0) < 1.0e-10 );
    return l0;
  }
  
public:
  std::vector<CRender2Tex_DrawOldGL> aSampler;
};

}
}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/render2tex_glold.cpp"
#endif

#endif /* depth_hpp */
