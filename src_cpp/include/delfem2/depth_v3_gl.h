#ifndef DEPTH_H
#define DEPTH_H

#include <stdio.h>

class CDepthContext
{
public:
  CDepthContext(){
    id_framebuffer = -1;
    id_depth_texture = -1;
  }
  CDepthContext(const std::vector<int>& winSize){
    this->Init(winSize[0],winSize[1]);
  }
  void Init(int width, int height){
    this->width = width;
    this->height = height;
    this->SetFrameBufferSize(width,height);
  }
  void DeleteFrameBuffer();
  void SetFrameBufferSize(int width, int height);
  void Start() const{
    glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, id_depth_render_buffer);
  }
  void End() const {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }
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
    color.resize(4);
    color[0] = 1;  color[1] = 0;  color[2] = 0;  color[3] = 1;
    nResW = 256;
    nResH = 256;
    lengrid = 0.01;
    orgPrj[  0]=0; orgPrj[  1]=0; orgPrj[  2]=0;
    dirPrj[  0]=0; dirPrj[  1]=0; dirPrj[  2]=1;
    dirWidth[0]=1; dirWidth[1]=0; dirWidth[2]=0;
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
  void Draw_Point() const;
  void SetView();
  void getGPos(double p[3], int i, double depth) const;
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
  int nResW;
  int nResH;
  double lengrid;
  double depth_max;
  double dirPrj[3];
  double dirWidth[3];
  double orgPrj[3];
  std::vector<float> aDepth;
  std::vector<double> color;
private:
  GLint view[4];
};

#endif /* depth_hpp */
