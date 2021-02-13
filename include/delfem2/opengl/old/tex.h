#ifndef DFM2_OPENGL_OLD_TEX_H
#define DFM2_OPENGL_OLD_TEX_H

#include "delfem2/opengl/tex.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/mat3.h"

namespace delfem2 {
namespace opengl {

void DrawTexMaxSize(
    const delfem2::opengl::CTexRGB &tex0) {
  ::glEnable(GL_TEXTURE_2D);
  ::glBindTexture(GL_TEXTURE_2D, tex0.id_tex);
  glColor3f(1.f, 1.f, 1.f);
  glBegin(GL_QUADS);
  ::glTexCoord2d(0, 0);
  glVertex2f(-1, +1);
  ::glTexCoord2d(1, 0);
  glVertex2f(+1, +1);
  ::glTexCoord2d(1, 1);
  glVertex2f(+1, -1);
  ::glTexCoord2d(0, 1);
  glVertex2f(-1, -1);
  glEnd();
}


void DrawCamInteriorMatrixWithTexture(
    const delfem2::opengl::CTexRGB &tex0,
    const float Kinv[9],
    float tex_loc_z) {
  ::glEnable(GL_TEXTURE_2D);
  const float ax[4][3] = {
      {0.f,           0.f,           1.f},
      {float(tex0.w), 0.f,           1.f},
      {float(tex0.w), float(tex0.h), 1.f},
      {0.f,           float(tex0.h), 1.f}};
  float aX[4][3];
  for (unsigned int i = 0; i < 4; ++i) {
    delfem2::MatVec3(aX[i], Kinv, ax[i]);
    const float a = tex_loc_z / aX[i][2];
    aX[i][0] *= a;
    aX[i][1] *= a;
    aX[i][2] = tex_loc_z;

  }
  ::glBindTexture(GL_TEXTURE_2D, tex0.id_tex);
  glColor3f(1.f, 1.f, 1.f);
  glBegin(GL_QUADS);
  ::glTexCoord2d(0, 0);
  glVertex3fv(aX[0]);
  ::glTexCoord2d(1, 0);
  glVertex3fv(aX[1]);
  ::glTexCoord2d(1, 1);
  glVertex3fv(aX[2]);
  ::glTexCoord2d(0, 1);
  glVertex3fv(aX[3]);
  glEnd();
  //
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_LINES);
  glVertex3fv(aX[0]);
  glVertex3d(0, 0, 0);
  glVertex3fv(aX[1]);
  glVertex3d(0, 0, 0);
  glVertex3fv(aX[2]);
  glVertex3d(0, 0, 0);
  glVertex3fv(aX[3]);
  glVertex3d(0, 0, 0);
  ::glEnd();
}

void DrawCam(
    const delfem2::opengl::CTexRGB &tex0,
    const float R[9],
    const float t[3],
    const float K[9]) {
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  float A[16];
  delfem2::Mat4_AffineTrans_RotTransl(
      A,
      R, t);
  ::glMultMatrixf(A);
  float Kinv[9];
  delfem2::Inverse_Mat3(Kinv, K);
  DrawCamInteriorMatrixWithTexture(tex0, Kinv, -1.f);
  delfem2::opengl::DrawAxis(1);
  ::glPopMatrix();
}

void DrawTex_CvNormDist(
    const std::vector<float> &aNormDist,
    const delfem2::opengl::CTexRGB &tex0,
    const float *R,
    const float *t,
    const float *K)
{
  const unsigned int nw = tex0.w;
  const unsigned int nh = tex0.h;
  assert(aNormDist.size()==nw*nh*4);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  float A[16];
  delfem2::Mat4_AffineTrans_RotTransl(
      A,
      R, t);
  ::glMultMatrixf(A);
  //
  double Kinv[9];
  delfem2::Inverse_Mat3(Kinv, K);
  double ratio0 = 1. / Kinv[8];

  float B[16];
  delfem2::Mat4_Mat3(B, Kinv);
  float C[16];
  delfem2::Transpose_Mat4(C, B);
  ::glMultMatrixf(C);
  //
  ::glDisable(GL_POLYGON_OFFSET_FILL);
  ::glEnable(GL_TEXTURE_2D);
  ::glBindTexture(GL_TEXTURE_2D, tex0.id_tex);
  ::glColor3d(1, 1, 1);
  ::glBegin(GL_QUADS);
  for (unsigned int iw = 0; iw < nw - 1; ++iw) {
    for (unsigned int ih = 0; ih < nh - 1; ++ih) {
      const float x0[3] = {float(iw)+0.5f,float(ih)+0.5f,1.f};
      const float x1[3] = {float(iw)+1.5f,float(ih)+0.5f,1.f};
      const float x2[3] = {float(iw)+1.5f,float(ih)+1.5f,1.f};
      const float x3[3] = {float(iw)+0.5f,float(ih)+1.5f,1.f};
      float Kix0[3]; MatVec3(Kix0, Kinv, x0);
      float Kix1[3]; MatVec3(Kix1, Kinv, x1);
      float Kix2[3]; MatVec3(Kix2, Kinv, x2);
      float Kix3[3]; MatVec3(Kix3, Kinv, x3);
      float ntKix0 = Dot3(aNormDist.data()+((ih + 0) * nw + (iw + 0))*4, Kix0);
      float ntKix1 = Dot3(aNormDist.data()+((ih + 0) * nw + (iw + 1))*4, Kix1);
      float ntKix2 = Dot3(aNormDist.data()+((ih + 1) * nw + (iw + 1))*4, Kix2);
      float ntKix3 = Dot3(aNormDist.data()+((ih + 1) * nw + (iw + 0))*4, Kix3);
//      std::cout << ntKix0 << " " << ntKix1 << " " << ntKix2 << " " << ntKix3 << std::endl;
      double z0 = -aNormDist[((ih + 0) * nw + (iw + 0))*4+3] / ntKix0;
      double z1 = -aNormDist[((ih + 0) * nw + (iw + 1))*4+3] / ntKix1;
      double z2 = -aNormDist[((ih + 1) * nw + (iw + 1))*4+3] / ntKix2;
      double z3 = -aNormDist[((ih + 1) * nw + (iw + 0))*4+3] / ntKix3;
      ::glTexCoord2d((iw + 0.5) / nw, (ih + 0.5) / nh);
      ::glVertex3d(z0 * (iw + 0.5), z0 * (ih + 0.5), z0);
      ::glTexCoord2d((iw + 1.5) / nw, (ih + 0.5) / nh);
      ::glVertex3d(z1 * (iw + 1.5), z1 * (ih + 0.5), z1);
      ::glTexCoord2d((iw + 1.5) / nw, (ih + 1.5) / nh);
      ::glVertex3d(z2 * (iw + 1.5), z2 * (ih + 1.5), z2);
      ::glTexCoord2d((iw + 0.5) / nw, (ih + 1.5) / nh);
      ::glVertex3d(z3 * (iw + 0.5), z3 * (ih + 1.5), z3);
    }
  }
  ::glEnd();
  //
  ::glPopMatrix();
}

}
}

#endif