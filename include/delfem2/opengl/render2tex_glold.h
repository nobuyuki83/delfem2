/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_GLOLD_SMPLR_H
#define DFM2_GLOLD_SMPLR_H

#include <stdio.h>
#include <vector>
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
  std::vector<CRender2Tex_DrawOldGL> aSampler;
  
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
  
};

}
}

#endif /* depth_hpp */
