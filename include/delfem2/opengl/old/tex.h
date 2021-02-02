#ifndef OPENGL_OLD_TEX_H
#define OPENGL_OLD_TEX_H

#include "delfem2/opengl/tex.h"

void DrawTexMaxSize(
    const delfem2::opengl::CTexRGB& tex0)
{
  ::glEnable(GL_TEXTURE_2D);
  ::glBindTexture(GL_TEXTURE_2D,tex0.id_tex);
  glColor3f(1.f, 1.f, 1.f);
  glBegin(GL_QUADS);
  ::glTexCoord2d(0,0); glVertex2f(-1,+1);
  ::glTexCoord2d(1,0); glVertex2f(+1,+1);
  ::glTexCoord2d(1,1); glVertex2f(+1,-1);
  ::glTexCoord2d(0,1); glVertex2f(-1,-1);
  glEnd();
}

#endif