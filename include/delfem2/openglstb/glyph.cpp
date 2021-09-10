
#include "delfem2/openglstb/glyph.h"

#include <iostream>
#include <fstream>

#if defined(_WIN32) // windows
  #include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
  #define GL_SILENCE_DEPRECATION
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#include "delfem2/str.h"

// --------------------------------------

 delfem2::openglstb::CGlyph::CGlyph(const std::string &fpath) {
  int channels;
  unsigned char *img = stbi_load(
      fpath.c_str(),
      &width, &height, &channels, 0);
  aRGBA.assign(img, img + width * height * channels);
}

void delfem2::openglstb::CGlyph::InitGL() {
  glEnable(GL_TEXTURE_2D);
  glGenTextures(1, &texName);
  glBindTexture(GL_TEXTURE_2D, texName);
  glTexImage2D(
      GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0,
      GL_RGBA, GL_UNSIGNED_BYTE, aRGBA.data());
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}


void delfem2::openglstb::CGlyph::DrawEntireGlyph() const {
  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  ::glBindTexture(GL_TEXTURE_2D, texName);
  ::glDisable(GL_LIGHTING);
  ::glEnable(GL_TEXTURE_2D);
  const double asp = width / double(height);
  //
  ::glBegin(GL_QUADS);
  ::glTexCoord2d(0, 1);
  ::glVertex3d(-asp, -1, 0);
  ::glTexCoord2d(1, 1);
  ::glVertex3d(+asp, -1, 0);
  ::glTexCoord2d(1, 0);
  ::glVertex3d(+asp, +1, 0);
  ::glTexCoord2d(0, 0);
  ::glVertex3d(-asp, +1, 0);
  ::glEnd();
}


void delfem2::openglstb::CGlyph::ParseGlyphInfo(const std::string &fpath) {
  std::ifstream fin(fpath);
  std::string line;
  getline(fin, line);
  getline(fin, line);
  getline(fin, line);
  getline(fin, line);
  auto aToken = Split(line, " =");
  assert(aToken.size() == 3);
  unsigned int nchar = std::strtol(aToken[2].data(), 0, 10);
  for (unsigned int ichar = 0; ichar < nchar; ++ichar) {
    getline(fin, line);
    auto aVal = Split(line, ' ');
    std::cout << line << std::endl;
    SData data;
	int id = 0;
    for (const std::string &s: aVal) {
      if (s.find("id=") == 0) { id = std::strtol(s.data() + 3, 0, 10); }
      if (s.find("x=") == 0) { data.x = std::strtol(s.data() + 2, 0, 10); }
      if (s.find("y=") == 0) { data.y = std::strtol(s.data() + 2, 0, 10); }
      if (s.find("width=") == 0) { data.width = std::strtol(s.data() + 6, 0, 10); }
      if (s.find("height=") == 0) { data.height = std::strtol(s.data() + 7, 0, 10); }
      if (s.find("xadvance=") == 0) { data.xadvance = std::strtol(s.data() + 9, 0, 10); }
    }
    mapData.insert(std::make_pair(
		static_cast<char>(id), 
		data));
  }
}

double delfem2::openglstb::CGlyph::DrawCharAt(char c, double scale, double px, double py) {
  auto itr = mapData.find(c);
  if (itr == mapData.end()) { return px; }
  const SData &data = itr->second;
  const double x0 = data.x / double(width);
  const double y0 = data.y / double(height);
  const double x1 = (data.x + data.width) / double(width);
  const double y1 = (data.y + data.height) / double(height);
  const double w0 = scale * data.width;
  const double h0 = scale * data.height;
  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  ::glBindTexture(GL_TEXTURE_2D, texName);
  ::glDisable(GL_LIGHTING);
  ::glEnable(GL_TEXTURE_2D);
  //
  ::glBegin(GL_QUADS);
  ::glTexCoord2d(x0, y1);
  ::glVertex2d(px, py);
  ::glTexCoord2d(x1, y1);
  ::glVertex2d(px + w0, py);
  ::glTexCoord2d(x1, y0);
  ::glVertex2d(px + w0, py + h0);
  ::glTexCoord2d(x0, y0);
  ::glVertex2d(px, py + h0);
  ::glEnd();
  return px + data.xadvance * scale;
}

void delfem2::openglstb::CGlyph::DrawStringAt(const std::string &str, double scale, double px, double py) {
  for (char c: str) {
    px = DrawCharAt(c, scale, px, py);
  }
}

#if defined(__APPLE__) && defined(__MACH__)
  #undef GL_SILENCE_DEPRECATION
#endif
