#ifndef DEPTH_H
#define DEPTH_H

#include <stdio.h>

class CInputDepth
{
public:
  virtual void Draw() const = 0;
};

class CGPUSampler
{
public:
  CGPUSampler(){
    nResX = 0;
    nResY = 0;
    isColor = false;
    isDepth = false;
    id_tex_color = 0;
    //////
    color.resize(4);
    color[0] = 1;  color[1] = 0;  color[2] = 0;  color[3] = 1;
    lengrid = 0.01;
    origin[  0]=0; origin[  1]=0; origin[  2]=0;
    z_axis[  0]=0; z_axis[  1]=0; z_axis[  2]=1;
    x_axis[0]=1; x_axis[1]=0; x_axis[2]=0;
    draw_len_axis = 1.0;
  }
  CGPUSampler(int nw, int nh, bool isColor, bool isDepth){
    this->Init(nw,nh,isColor,isDepth);
  }
  void Init(int nw, int nh, bool isColor, bool isDepth);
  void LoadTex();
  /////
  void Draw() const;
  std::vector<double> MinMaxXYZ() const {
    std::vector<double> mm(6);
    mm[0] = +1;
    mm[1] = -1;
    return mm;
  }
  /////
  void Draw_Axis() const;
  void Draw_Point() const;
  void Draw_BoundingBox() const;
  void SetView();
  void getGPos(double p[3], int ix, int iy, double depth) const;
  ////
  void SetColor(double r, double g, double b);
  void SaveDepthCSV(const std::string& path) const;
  void SetCoord(double elen, double depth_max,
                const std::vector<double>& orgPrj,
                const std::vector<double>& dirPrj,
                const std::vector<double>& dirWidth);
  void Start();
  void End();
public:
  bool isColor;
  bool isDepth;
  int nResX;
  int nResY;
  double lengrid;
  double z_range;
  double z_axis[3];
  double x_axis[3];
  double origin[3];
  std::vector<float> aZ;
  std::vector<unsigned char> aRGBA;
  /////
  std::vector<double> color;
  double draw_len_axis;
  unsigned int id_tex_color;
private:
  GLint view[4];
};

#endif /* depth_hpp */
