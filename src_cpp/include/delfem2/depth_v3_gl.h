#ifndef DEPTH_H
#define DEPTH_H

#include <stdio.h>

class CInputDepth
{
public:
  virtual void Draw() const = 0;
};

class CDepth
{
public:
  CDepth(){
    color.resize(4);
    color[0] = 1;  color[1] = 0;  color[2] = 0;  color[3] = 1;
    nResX = 256;
    nResY = 256;
    lengrid = 0.01;
    origin[  0]=0; origin[  1]=0; origin[  2]=0;
    z_axis[  0]=0; z_axis[  1]=0; z_axis[  2]=1;
    x_axis[0]=1; x_axis[1]=0; x_axis[2]=0;
    draw_len_axis = 1.0;
  }
  CDepth(int nw, int nh, double l, double dm,
         const std::vector<double>& org,
         const std::vector<double>& dirP,
         const std::vector<double>& dirW){
    this->SetCoord(nw,nh,l,dm,org,dirP,dirW);
  }
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
  void SetCoord(int nresw, int nresh, double elen,
                double depth_max,
                const std::vector<double>& orgPrj,
                const std::vector<double>& dirPrj,
                const std::vector<double>& dirWidth);
  void Start();
  void End();
public:
  int nResX;
  int nResY;
  double lengrid;
  double z_range;
  double z_axis[3];
  double x_axis[3];
  double origin[3];
  std::vector<float> aZ;
  /////
  std::vector<double> color;
  double draw_len_axis;
private:
  GLint view[4];
};

#endif /* depth_hpp */
