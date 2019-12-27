/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef tex_h
#define tex_h

#include <vector>
#include <string>

class CTexture
{
public:
  std::vector<unsigned char> aRGB;
  unsigned int id_tex;
  unsigned int h,w;
  double min_x, max_x, min_y, max_y;
  
public:
  CTexture(){
    id_tex = 0;
    this->w = 0;
    this->h = 0;
    aRGB.clear();
  }
  CTexture(int w, int h, const unsigned char* pD, const std::string& typeData){
    this->h = h;
    this->w = w;
    this->aRGB.assign(pD,pD+h*w*3);
    if( typeData == "bgr" ){
      for(int i=0;i<h*w;++i){ // rgb -> bgr
        unsigned char b0 = aRGB[i*3+0];
        unsigned char r0 = aRGB[i*3+2];
        aRGB[i*3+0] = r0;
        aRGB[i*3+2] = b0;
      }
    }
    id_tex = 0;
    this->min_x = 0.0;
    this->max_x = (double)w;
    this->min_y = 0.0;
    this->max_y = (double)h;
  }
  
  void LoadTex();
  
  void Draw();
  
  std::vector<double>  MinMaxXYZ() const {
    std::vector<double> m(6,0.0);
    m[0] = this->min_x;
    m[1] = this->max_x;
    m[2] = this->min_y;
    m[3] = this->max_y;
    m[4] = 0;
    m[5] = 0;
    return m;
  }
  void SetMinMaxXY(const std::vector<double>& mmxy){
    if( mmxy.size() < 4 ){ return; }
    this->min_x = mmxy[0];
    this->max_x = mmxy[1];
    this->min_y = mmxy[2];
    this->max_y = mmxy[3];
  }
};


// ---------------------------------------
// texture related funcitons

void LoadImage_PPM(const std::string& filename,
                   std::vector<unsigned char>& image,
                   int& width, int& height);

void SaveImage(const std::string& path);


class SFile_TGA
{
public:
  unsigned char imageTypeCode;
  short int imageWidth;
  short int imageHeight;
  unsigned char bitCount;
  unsigned char *imageData;
};

bool LoadTGAFile(const std::string& filename, SFile_TGA *tgaFile);
int ReadPPM_SetTexture(const std::string& fname);

unsigned int LoadTexture(
    const unsigned char* image,
    const int width, const int height, const int bpp);
void DrawTextureBackground(
    const unsigned int tex,
     const int imgWidth,
     const int imgHeight,
     const int winWidth,
     const int winHeight);

class CTextureInfo
{
public:
  std::string full_path;
  int width, height, bpp; // byte par pixel
  int id_tex_gl;
};

class CTexManager
{
public:
  void Clear();
  void AddTexture(const unsigned char* pixels,
                  const std::string& path,
                  int width, int height, int bpp)
  {
    const int id_tex_gl = LoadTexture(pixels, width,height,bpp);
    CTextureInfo texinfo;
    texinfo.full_path = path;
    texinfo.height = height;
    texinfo.width = width;
    texinfo.bpp = bpp;
    texinfo.id_tex_gl = id_tex_gl;
    /////
    bool is_new = true;
    for(int itex=0;itex<(int)aTexInfo.size();++itex){
      if( aTexInfo[itex].full_path != path ) continue;
      aTexInfo[itex] = texinfo;
      is_new = false;
    }
    if( is_new ){
      aTexInfo.push_back(texinfo);
    }
  }
  void AddPath(const std::string& path){
    CTextureInfo tex;
    tex.width = -1;
    tex.height = -1;
    tex.bpp = -1;
    tex.id_tex_gl = -1;
    tex.full_path = path;
    aTexInfo.push_back(tex);
  }
  void BindTexturePath(const std::string& path) const;
public:
  std::vector<CTextureInfo> aTexInfo;
};

#endif
