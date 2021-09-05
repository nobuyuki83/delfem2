/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/tex.h"
#include <fstream>
#include <iostream>
#include <cassert>
#include <cstdlib>

#if defined(_MSC_VER)
#  pragma warning( push )
// C4996 (lev3): Your code uses a function, class member, variable, or typedef that's marked deprecated.
#  pragma warning( disable : 4996 )
#elif defined(__GNUC__) || defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#if defined(_WIN32) // windows
#  include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__) // mac
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif


// ------------------------

void delfem2::opengl::SaveViewportAsImagePpm(
    const std::string &path) {
  static unsigned int inum = 0;
  int viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  void *image = malloc(3 * viewport[2] * viewport[3]);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(
      0, 0,
      viewport[2], viewport[3],
      GL_RGB, GL_UNSIGNED_BYTE, image);
  unsigned int width = viewport[2];
  unsigned int height = viewport[3];
  std::ofstream fout;
  //  //  fname << "out";
  fout.open(path.c_str(), std::ios::out);
  fout << "P3\n";
  fout << width << " " << height << "\n";
  fout << "255\n";
  //  fout << "255\n";
  //  fout << "255\n";
  char *img = (char *) image;
  for (unsigned int ih = 0; ih < height; ih++) {
    for (unsigned int iw = 0; iw < width; iw++) {
      unsigned int i = (height - 1 - ih) * width + iw;
      int r = (unsigned char) img[i * 3 + 0];
      int g = (unsigned char) img[i * 3 + 1];
      int b = (unsigned char) img[i * 3 + 2];
      fout << r << " " << g << " " << b << "\n";
      //    std::cout << i << " " << r << " "<< g << " "<< b << std::endl;
    }
  }
  fout.close();
  //  if( inum >= 700 ) abort();
  //  if( inum >= 400 ) abort();
  if (inum >= 600) abort();
  inum++;
}

int delfem2::opengl::SetTexture_RGB(
    unsigned int w, unsigned int h,
    const std::vector<unsigned char> &image) {
  glEnable(GL_TEXTURE_2D);
  GLuint m_texName = 0;
  glGenTextures(1, &m_texName);
  glBindTexture(GL_TEXTURE_2D, m_texName);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
               static_cast<int>(w),
               static_cast<int>(h),
               0, GL_RGB, GL_UNSIGNED_BYTE, image.data());

  return (int) m_texName;
}

GLuint delfem2::opengl::LoadTexture
    (const unsigned char *image,
     const int width, const int height, const int bpp) {
  GLuint id_tex;
  glGenTextures(1, &id_tex);
  glBindTexture(GL_TEXTURE_2D, id_tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  if (bpp == 3) {
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D,
                 0, GL_RGBA, (GLsizei) width, (GLsizei) height,
                 0, GL_RGB, GL_UNSIGNED_BYTE, image);
  }
  if (bpp == 4) {
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glTexImage2D(GL_TEXTURE_2D,
                 0, GL_RGBA, (GLsizei) width, (GLsizei) height,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
  }
  return id_tex;
};

// use it for GLSL shader drawing
void DrawRectangle_FullCanvas_oldGL() {
#ifdef GL_MODELVIEW
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glBegin(GL_QUADS);
  glVertex2d(-1, -1);
  glVertex2d(+1, -1);
  glVertex2d(+1, +1);
  glVertex2d(-1, +1);
  glEnd();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glPopAttrib();
#endif
}

// ================================================================

void delfem2::opengl::CTexRGB::InitGL() {
  if (id_tex == 0) { ::glGenTextures(1, &id_tex); }
  glBindTexture(GL_TEXTURE_2D, id_tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  assert(pixel_color.size() == width * height * channels );
  if( this->channels ==  3 ) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
                 static_cast<int>(width),
                 static_cast<int>(height),
                 0, GL_RGB, GL_UNSIGNED_BYTE,
                 pixel_color.data());
  }
  else if( this->channels ==  4 ) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                 static_cast<int>(width),
                 static_cast<int>(height),
                 0, GL_RGBA, GL_UNSIGNED_BYTE,
                 pixel_color.data());
  }
  glBindTexture(GL_TEXTURE_2D, 0);
}





// --------------------------------------------------


void delfem2::opengl::CTexRGB_Rect2D::Draw_oldGL() const {
  if (id_tex == 0) { return; }
#ifdef GL_LIGHTING
  ::glEnable(GL_TEXTURE_2D);
  ::glDisable(GL_LIGHTING);
  ::glBindTexture(GL_TEXTURE_2D, id_tex);
  ::glColor3d(1, 1, 1);
  ::glBegin(GL_QUADS);
  ::glTexCoord2d(0.0, 0.0);
  ::glVertex3d(min_x, min_y, z);
  ::glTexCoord2d(1.0, 0.0);
  ::glVertex3d(max_x, min_y, z);
  ::glTexCoord2d(1.0, 1.0);
  ::glVertex3d(max_x, max_y, z);
  ::glTexCoord2d(0.0, 1.0);
  ::glVertex3d(min_x, max_y, z);
  ::glEnd();
  ::glBindTexture(GL_TEXTURE_2D, 0);
  ::glDisable(GL_TEXTURE_2D);
#endif
}

/*
void DrawTextureBackground
(const GLuint tex,
 const int imgWidth,
 const int imgHeight,
 const int winWidth,
 const int winHeight)
{
  double imgAsp = (double)imgWidth/imgHeight;
  double winAsp = (double)winWidth/winHeight;
  /////
  glPushAttrib(GL_TRANSFORM_BIT|GL_CURRENT_BIT|GL_ENABLE_BIT);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  
  if (winAsp>imgAsp){
    double imgWidth2 = imgHeight*winAsp;
    gluOrtho2D(
               -0.5*imgWidth2,+0.5*imgWidth2,
               -0.5*imgHeight,+0.5*imgHeight);
  }
  else{
    double imgHeight2 = (double)imgWidth/winAsp;
    gluOrtho2D(
               -0.5*imgWidth, +0.5*imgWidth,
               -0.5*imgHeight2, +0.5*imgHeight2);
  }
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  
  glBindTexture(GL_TEXTURE_2D, tex);
  
  glColor3f(1.0f, 1.0f, 1.0f);
  glBegin(GL_QUADS);
  glTexCoord2i(0, 0); glVertex2d(-0.5*imgWidth, -0.5*imgHeight);
  glTexCoord2i(1, 0); glVertex2d(+0.5*imgWidth, -0.5*imgHeight);
  glTexCoord2i(1, 1); glVertex2d(+0.5*imgWidth, +0.5*imgHeight);
  glTexCoord2i(0, 1); glVertex2d(-0.5*imgWidth, +0.5*imgHeight);
  glEnd();
  
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  
  glPopAttrib();
}
 */

void delfem2::opengl::CTexManager::Clear() {
  for (auto &itex : aTexInfo) {
    unsigned int id_tex_gl = itex.id_tex_gl;
    if (glIsTexture(id_tex_gl)) {
      ::glDeleteTextures(1, &id_tex_gl);
    }
  }
  aTexInfo.clear();
}

void delfem2::opengl::CTexManager::BindTexturePath(const std::string &path) const {
  for (const auto &itex : aTexInfo) {
    if (itex.full_path != path) continue;
    glBindTexture(GL_TEXTURE_2D, itex.id_tex_gl);
    glEnable(GL_TEXTURE_2D);
  }
}

#if defined(_MSC_VER)
#pragma warning( pop )
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
