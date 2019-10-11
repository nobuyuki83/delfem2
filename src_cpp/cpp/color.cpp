/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstring>
#include <cstdlib>

#include "delfem2/color.h"

void GetRGB_HSV
(float&r, float& g, float& b,
 float h, float s, float v)
{
  r = v;
  g = v;
  b = v;
  if (s > 0.0f) {
    h *= 6.0f;
    const int i = (int) h;
    const float f = h - (float) i;
    switch (i) {
      default:
      case 0:
        g *= 1 - s * (1 - f);
        b *= 1 - s;
        break;
      case 1:
        r *= 1 - s * f;
        b *= 1 - s;
        break;
      case 2:
        r *= 1 - s;
        b *= 1 - s * (1 - f);
        break;
      case 3:
        r *= 1 - s;
        g *= 1 - s * f;
        break;
      case 4:
        r *= 1 - s * (1 - f);
        g *= 1 - s;
        break;
      case 5:
        g *= 1 - s;
        b *= 1 - s * f;
        break;
    }
  }
}

// ------------------------------------------------------------

void interpolateColor
(CColor& Cout, float r, const CColor& C0, const CColor& C1)
{
  Cout.r = (1-r)*C0.r+r*C1.r;
  Cout.g = (1-r)*C0.g+r*C1.g;
  Cout.b = (1-r)*C0.b+r*C1.b;
  Cout.a = (1-r)*C0.a+r*C1.a;
}

void heatmap(double input,double* color)
{
  if(0 <=input&&input <=0.25){
    color[0] = 0.0;
    color[1] = input*4.0;
    color[2] = 1.0;
  }
  else if(0.25<input && input <=0.5){
    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 2.0-input*4.0;
  }
  else if(0.5<input && input <=0.75){
    color[0] = input*4 -2.0;
    color[1] = 1.0;
    color[2] = 0.0;
  }
  else if(0.75<input&& input <=1.0){
    color[0] = 1.0;
    color[1] = 4.0 -input*4.0;
    color[2] = 0.0;
  }
  else if(1.0<input){
    color[0] = 1.0;
    color[1] = 0.0;
    color[2] = 0.0;
  }
  else{
    color[0] = 0.0;
    color[1] = 0.0;
    color[2] = 1.0;
  }
}

CColor getColor(double input, const std::vector<std::pair<double, CColor> >& colorMap)
{
  if (colorMap.size()==0) return CColor::Black();
  if (input < colorMap[0].first){
    return colorMap[0].second;
  }
  for (int ic = 0; ic<(int)colorMap.size()-1; ++ic){
    double val0 = colorMap[ic].first;
    double val1 = colorMap[ic+1].first;
    if (val0<=input&&input<=val1){
      float rp = (float)((input-val0)/(val1-val0));
      CColor color;
      interpolateColor(color, rp, colorMap[ic].second, colorMap[ic+1].second);
      return color;
    }
  }
  return colorMap[colorMap.size()-1].second;
}

void makeHeatMap_BlueGrayRed(std::vector<std::pair<double, CColor> >& colorMap, float min, float max)
{
  double diff = (max-min)*0.25;
  colorMap.push_back(std::make_pair(min+diff*0, CColor(0.0f, 0.0f, 1.0f, 1.0f))); // blue
  colorMap.push_back(std::make_pair(min+diff*1, CColor(0.0f, 0.2f, 1.0f, 1.0f)));
  colorMap.push_back(std::make_pair(min+diff*2, CColor(0.5f, 0.5f, 0.5f, 1.0f))); // gray
  colorMap.push_back(std::make_pair(min+diff*3, CColor(1.0f, 0.2f, 0.0f, 1.0f)));
  colorMap.push_back(std::make_pair(min+diff*4, CColor(1.0f, 0.0f, 0.0f, 1.0f))); // red
}

void makeHeatMap_BlueCyanGreenYellowRed(std::vector<std::pair<double, CColor> >& colorMap, float min, float max, float alpha)
{
  double diff = (max-min)*0.25;
  colorMap.push_back(std::make_pair(min+diff*0, CColor(0.0f, 0.0f, 1.0f, alpha))); // blue
  colorMap.push_back(std::make_pair(min+diff*1, CColor(0.0f, 1.0f, 1.0f, alpha))); // cyan
  colorMap.push_back(std::make_pair(min+diff*2, CColor(0.0f, 1.0f, 0.0f, alpha))); // green
  colorMap.push_back(std::make_pair(min+diff*3, CColor(1.0f, 1.0f, 0.0f, alpha))); // yellow
  colorMap.push_back(std::make_pair(min+diff*4, CColor(1.0f, 0.0f, 0.0f, alpha))); // red
}

void makeHeatMap_RedYellowGreenCyanBlue(std::vector<std::pair<double, CColor> >& colorMap, float min, float max)
{
  double diff = (max-min)*0.25;
  colorMap.push_back(std::make_pair(min+diff*0, CColor(1.0f, 0.0f, 0.0f, 1.0f))); // red
  colorMap.push_back(std::make_pair(min+diff*1, CColor(1.0f, 1.0f, 0.0f, 1.0f))); // yellow
  colorMap.push_back(std::make_pair(min+diff*2, CColor(0.0f, 1.0f, 0.0f, 1.0f))); // green
  colorMap.push_back(std::make_pair(min+diff*3, CColor(0.0f, 1.0f, 1.0f, 1.0f))); // cyan
  colorMap.push_back(std::make_pair(min+diff*4, CColor(0.0f, 0.0f, 1.0f, 1.0f))); // blue
}

// ----------------------------------------------------

void Write_Ply_Tri2DMesh_HeightColor
(const std::string& fname,
 const std::vector<int>& aTri1,
 const std::vector<double>& aXY1,
 const std::vector<double>& aVal,
 std::vector< std::pair<double,CColor> >& colorMap)
{
  const int np = aXY1.size()/2;
  const int ntri = aTri1.size()/3;
  std::ofstream fout;
  fout.open(fname.c_str(),std::ios::out);
  fout << "ply" << std::endl;
  fout << "format ascii 1.0" << std::endl;
  fout << "element vertex " << np << std::endl;
  fout << "property float x" << std::endl;
  fout << "property float y" << std::endl;
  fout << "property float z" << std::endl;
  fout << "property uchar red" << std::endl;
  fout << "property uchar green" << std::endl;
  fout << "property uchar blue" << std::endl;
  fout << "element face " << ntri << std::endl;
  fout << "property list uchar int vertex_indices" << std::endl;
  fout << "end_header" << std::endl;
  for(int ip=0;ip<np;++ip){
    double v = aVal[ip];
    CColor c = getColor(v,colorMap);
    int cr, cg, cb;
    c.getRGBChar(cr,cg,cb);
    fout << aXY1[ip*2+0] << " " << aXY1[ip*2+1] << " " << v << " ";
    fout << cr << " " << cg << " " << cb << std::endl;
  }
  for(int itri=0;itri<ntri;++itri){
    fout << "3 " << aTri1[itri*3+0] << " " << aTri1[itri*3+1] << " " << aTri1[itri*3+2] << std::endl;
  }
}


