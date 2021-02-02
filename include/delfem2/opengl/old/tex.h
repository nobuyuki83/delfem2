#ifndef OPENGL_OLD_TEX_H
#define OPENGL_OLD_TEX_H

#include "delfem2/opengl/tex.h"
#include "delfem2/opengl/old/funcs.h"

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


void DrawCamInteriorMatrixWithTexture(
    const delfem2::opengl::CTexRGB& tex0,
    const float Kinv[9],
    float tex_loc_z)
{
  ::glEnable(GL_TEXTURE_2D);
  const float ax[4][3] = {
      { 0.f,                     0.f, 1.f },
      { float(tex0.w),           0.f, 1.f },
      { float(tex0.w), float(tex0.h), 1.f },
      { 0.f,           float(tex0.h), 1.f } };
  float aX[4][3];
  for(unsigned int i=0;i<4;++i){
    delfem2::MatVec3(aX[i], Kinv, ax[i]);
    const float a = tex_loc_z / aX[i][2];
    aX[i][0] *= a;
    aX[i][1] *= a;
    aX[i][2] = tex_loc_z;

  }
  ::glBindTexture(GL_TEXTURE_2D,tex0.id_tex);
  glColor3f(1.f, 1.f, 1.f);
  glBegin(GL_QUADS);
  ::glTexCoord2d(0,0); glVertex3fv(aX[0]);
  ::glTexCoord2d(1,0); glVertex3fv(aX[1]);
  ::glTexCoord2d(1,1); glVertex3fv(aX[2]);
  ::glTexCoord2d(0,1); glVertex3fv(aX[3]);
  glEnd();
  //
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  glVertex3fv(aX[0]); glVertex3d(0, 0, 0);
  glVertex3fv(aX[1]); glVertex3d(0, 0, 0);
  glVertex3fv(aX[2]); glVertex3d(0, 0, 0);
  glVertex3fv(aX[3]); glVertex3d(0, 0, 0);
  ::glEnd();
}

void DrawCam(
    const delfem2::opengl::CTexRGB& tex0,
    const float R[9],
    const float t[3],
    const float K[9])
{
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  float A[16];
  delfem2::Mat4_AffineTrans_RotTransl(
      A,
      R,t);
  ::glMultMatrixf(A);
  float Kinv[9];
  delfem2::Inverse_Mat3(Kinv,K);
  DrawCamInteriorMatrixWithTexture(tex0, Kinv, -1.f);
  delfem2::opengl::DrawAxis(1);
  ::glPopMatrix();
}


#endif