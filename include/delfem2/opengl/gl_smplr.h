/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_GL_SMPLR_H
#define DFM2_GL_SMPLR_H

#include <stdio.h>
#include <vector>

class CGPUSampler
{
public:
  CGPUSampler(){
    nResX = 0;
    nResY = 0;
    id_tex_color = 0;
    bgcolor.resize(4);  bgcolor[0] = 1;  bgcolor[1] = 1;  bgcolor[2] = 1;  bgcolor[3] = 1;
    lengrid = 0.01;
    origin[  0]=0; origin[  1]=0; origin[  2]=0;
    z_axis[  0]=0; z_axis[  1]=0; z_axis[  2]=1;
    x_axis[0]=1; x_axis[1]=0; x_axis[2]=0;
  }
  CGPUSampler(int nw, int nh, std::string sFormatPixelColor, bool isDepth){
    this->Init(nw,nh,sFormatPixelColor,isDepth);
  }
  void InitGL();
  std::vector<double> MinMaxXYZ() const {
    std::vector<double> mm(6);
    mm[0] = +1;
    mm[1] = -1;
    return mm;
  }
  // ----------------------
  void Init(int nw, int nh, std::string sFormatPixelColor, bool isDepth);
  void Matrix_MVP(float mMV[16], float p[16]) const;
  std::vector<double> getGPos(int ix, int iy) const;
  void SaveDepthCSV(const std::string& path) const;
  void SetCoord(double elen, double depth_max,
                const std::vector<double>& orgPrj,
                const std::vector<double>& dirPrj,
                const std::vector<double>& dirWidth);
  void Start();
  void End();
  void ExtractFromTexture_Depth(std::vector<float>& aZ);
  void ExtractFromTexture_Color(std::vector<std::uint8_t>& aRGBA);
public:
  unsigned int id_tex_color;
  unsigned int id_tex_depth;
  unsigned int id_framebuffer;
  unsigned int nResX;
  unsigned int nResY;
  double lengrid;
  double z_range;
  double z_axis[3];
  double x_axis[3];
  double origin[3];
  std::vector<double> bgcolor;
  // --------------
protected:
  int view[4];
};

#endif /* depth_hpp */
