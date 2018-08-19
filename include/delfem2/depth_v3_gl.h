#ifndef DEPTH_H
#define DEPTH_H

#include <stdio.h>

#include "delfem2/vec3.h"

class CDepthContext
{
public:
  CDepthContext(){
    id_framebuffer = -1;
    id_depth_texture = -1;
  }
  void Init(int width, int height){
    this->width = width;
    this->height = height;
    this->SetFrameBufferSize(width,height);
  }
  void DeleteFrameBuffer();
  void SetFrameBufferSize(int width, int height);
public:
  unsigned int id_framebuffer;
  unsigned int id_depth_texture;
  unsigned int id_depth_render_buffer;
  int width;
  int height;
};

class CInputDepth
{
public:
  virtual void Draw() const = 0;
};

class CDepth
{
public:
  CDepth(){
    color[0] = 0;  color[1] = 0;  color[2] = 0;  color[3] = 1;
  }
  void Draw_Point(bool is_draw_miss) const;
  void TakeDepthShot(const CInputDepth& obj);
  void SetView();
  void getGPos(CVector3& p, int i, double depth) const;
  ////
  void SetColor(double r, double g, double b);
  void SaveDepthCSV(const std::string& path) const;
public:
  int nResW;
  int nResH;
  double lengrid;
  double depth_max;
  CVector3 dirPrj;
  CVector3 dirWidth;
  CVector3 orgPrj;
  std::vector<float> aDepth;
  double color[4];
};

#endif /* depth_hpp */
