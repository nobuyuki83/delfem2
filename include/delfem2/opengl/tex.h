/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * - this code should be compilable for both modern and legacy OpenGL.
 * - This code should not have any dependencies.
 */

#ifndef DFM2_OPENGL_TEX_H
#define DFM2_OPENGL_TEX_H

#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {
namespace opengl {

void SaveViewportAsImagePpm(const std::string &path);

int SetTexture_RGB(
    unsigned int w,
    unsigned int h,
    const std::vector<unsigned char> &image);

/**
 * @param bpp byte par pixel
 * @return texture id
 */
unsigned int LoadTexture(
    const unsigned char *image,
    int width, int height, int bpp);

// ----------------------------------

class CTexRGB {
 public:
  std::vector<unsigned char> aRGB;
  unsigned int id_tex_;
  unsigned int height_, width_;

 public:
  CTexRGB() {
    id_tex_ = 0;
    this->width_ = 0;
    this->height_ = 0;
    aRGB.clear();
  }

  /*
  CTexIdAndRGBData(int w, int h,
                   const unsigned char *pD,
                   const std::string &typeData) {
    this->Initialize(w,h, pD, typeData);
  }
   */

  virtual void Initialize(
      unsigned int w0,
      unsigned int h0,
      const unsigned char *pD,
      const std::string &typeData) {
    this->height_ = h0;
    this->width_ = w0;
    this->aRGB.assign(pD, pD + height_ * width_ * 3);
    if (typeData == "bgr") {
      for (unsigned int i = 0; i < height_ * width_; ++i) { // rgb -> bgr
        unsigned char b0 = aRGB[i * 3 + 0];
        unsigned char r0 = aRGB[i * 3 + 2];
        aRGB[i * 3 + 0] = r0;
        aRGB[i * 3 + 2] = b0;
      }
    }
    id_tex_ = 0;
  }

  virtual void InitGL();
};

// ----------------------------------

class CTexRGB_Rect2D :
    public CTexRGB {
 public:
  CTexRGB_Rect2D() : CTexRGB() {}
  virtual ~CTexRGB_Rect2D() = default;

  void Draw_oldGL() const;

  void Initialize(unsigned int w, unsigned int h,
                  const unsigned char *pD,
                  const std::string &typeData) override {
    CTexRGB::Initialize(w, h, pD, typeData);
    this->min_x = 0.0;
    this->max_x = (double) w;
    this->min_y = 0.0;
    this->max_y = (double) h;
  }

  std::vector<double> MinMaxAABB() const {
    return {this->min_x, this->min_y, z,
            this->max_x, this->max_y, z};
  }

  void SetMinMaxXY(const std::vector<double> &mmxy) {
    if (mmxy.size() < 4) { return; }
    this->min_x = mmxy[0];
    this->max_x = mmxy[1];
    this->min_y = mmxy[2];
    this->max_y = mmxy[3];
    z = (mmxy[4] + mmxy[5]) * 0.5;
  }
 public:
  // this is a coordinate for OpenGL image plane (after ModelView and Projection)
  double min_x = -1, max_x = +1;
  double min_y = -1, max_y = +1;
  double z = -1;
};

// ----------------------------------

class CTextureInfo {
 public:
  std::string full_path;
  int width = 0, height = 0, bpp = 0; // byte par pixel
  int id_tex_gl = 0;
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
#include "delfem2/opengl/tex.cpp"
#endif

#endif
