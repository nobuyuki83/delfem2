/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_TEX_GL_H
#define DFM2_TEX_GL_H

#include <vector>
#include <string>
#include "delfem2/dfm2_inline.h"

namespace delfem2 {
namespace opengl {

void SaveViewportAsImagePpm(const std::string &path);

int SetTexture_RGB
 (unsigned int w, unsigned int h,
  const std::vector<unsigned char> &image);

/**
 * @param bpp byte par pixel
 * @return texture id
 */
unsigned int LoadTexture
 (const unsigned char *image,
  const int width, const int height, const int bpp);

// ----------------------------------

class CTexRGB {
public:
  std::vector<unsigned char> aRGB;
  unsigned int id_tex;
  unsigned int h, w;
//  double min_x, max_x, min_y, max_y;

public:
  CTexRGB() {
    id_tex = 0;
    this->w = 0;
    this->h = 0;
    aRGB.clear();
  }

  /*
  CTexIdAndRGBData(int w, int h,
                   const unsigned char *pD,
                   const std::string &typeData) {
    this->Initialize(w,h, pD, typeData);
  }
   */
  
  virtual void Initialize(int w, int h,
                  const unsigned char *pD,
                  const std::string &typeData){
    this->h = h;
    this->w = w;
    this->aRGB.assign(pD, pD + h * w * 3);
    if (typeData == "bgr") {
      for (int i = 0; i < h * w; ++i) { // rgb -> bgr
        unsigned char b0 = aRGB[i * 3 + 0];
        unsigned char r0 = aRGB[i * 3 + 2];
        aRGB[i * 3 + 0] = r0;
        aRGB[i * 3 + 2] = b0;
      }
    }
    id_tex = 0;
  }

  virtual void InitGL();
};

// ----------------------------------

class CTexRGB_Rect2D :
public CTexRGB
{
public:
  CTexRGB_Rect2D() : CTexRGB() {
    min_x = min_y = 0;
    max_x = max_y = 0;
    z = 0;
  }
  
  void Draw_oldGL();
  
  virtual void Initialize(int w, int h,
                  const unsigned char *pD,
                  const std::string &typeData){
    CTexRGB::Initialize(w, h, pD, typeData);
    this->min_x = 0.0;
    this->max_x = (double) w;
    this->min_y = 0.0;
    this->max_y = (double) h;
  }
  
  std::vector<double> MinMaxXYZ() const {
    std::vector<double> m(6, 0.0);
    m[0] = this->min_x;
    m[1] = this->max_x;
    m[2] = this->min_y;
    m[3] = this->max_y;
    m[4] = z;
    m[5] = z;
    return m;
  }
  
  void SetMinMaxXY(const std::vector<double> &mmxy) {
    if (mmxy.size() < 4) { return; }
    this->min_x = mmxy[0];
    this->max_x = mmxy[1];
    this->min_y = mmxy[2];
    this->max_y = mmxy[3];
    z = (mmxy[4] + mmxy[5])*0.5;
  }
public:
  double min_x, max_x, min_y, max_y;
  double z;
};

// ----------------------------------

class CTextureInfo {
public:
  std::string full_path;
  int width, height, bpp; // byte par pixel
  int id_tex_gl;
};

class CTexManager {
public:
  void Clear();

  void AddTexture(const unsigned char *pixels,
                  const std::string &path,
                  int width, int height, int bpp) {
    const int id_tex_gl = LoadTexture(pixels, width, height, bpp);
    CTextureInfo texinfo;
    texinfo.full_path = path;
    texinfo.height = height;
    texinfo.width = width;
    texinfo.bpp = bpp;
    texinfo.id_tex_gl = id_tex_gl;
    /////
    bool is_new = true;
    for (int itex = 0; itex < (int) aTexInfo.size(); ++itex) {
      if (aTexInfo[itex].full_path != path) continue;
      aTexInfo[itex] = texinfo;
      is_new = false;
    }
    if (is_new) {
      aTexInfo.push_back(texinfo);
    }
  }

  void AddPath(const std::string &path) {
    CTextureInfo tex;
    tex.width = -1;
    tex.height = -1;
    tex.bpp = -1;
    tex.id_tex_gl = -1;
    tex.full_path = path;
    aTexInfo.push_back(tex);
  }

  void BindTexturePath(const std::string &path) const;

public:
  std::vector<CTextureInfo> aTexInfo;
};

}
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/tex_gl.cpp"
#endif

#endif
