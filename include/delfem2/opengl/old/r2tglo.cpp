/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#if defined(_WIN32) // windows
#  include <windows.h>
#endif
#include "glad/glad.h" // gl3.0+
#if defined(__APPLE__) && defined(__MACH__) // Mac
#  define GL_SILENCE_DEPRECATION
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/r2tglo.h"
#include "delfem2/opengl/funcs.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

// ------------------

DFM2_INLINE void delfem2::opengl::SetView(const CRender2Tex &r2t) {
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  ::glMultMatrixd(TransposeMat4ForOpenGL(r2t.mat_modelview,false).data());
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::glMultMatrixd(TransposeMat4ForOpenGL(r2t.mat_projection,true).data());
  ::glMatrixMode(GL_MODELVIEW);
}

// --------------------------------------------

void delfem2::opengl::CDrawerOldGL_Render2Tex::SetPointColor(double r, double g, double b) {
  colorPoint[0] = r;
  colorPoint[1] = g;
  colorPoint[2] = b;
}

void delfem2::opengl::CDrawerOldGL_Render2Tex::Draw_Texture(const CRender2Tex &r2t) {
  double mMVP[16];
  MatMat4(
      mMVP,
      r2t.mat_projection,
      r2t.mat_modelview);
  double mMVPinv[16];
  Inverse_Mat4(mMVPinv, mMVP);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glMultMatrixd(TransposeMat4ForOpenGL(mMVPinv,false).data());
  //
  ::glEnable(GL_TEXTURE_2D);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1, 1, 1);
  ::glBegin(GL_QUADS);
  ::glTexCoord2d(0.0, 0.0);
  ::glVertex3d(-1, -1, +1);
  ::glTexCoord2d(1.0, 0.0);
  ::glVertex3d(+1, -1, +1);
  ::glTexCoord2d(1.0, 1.0);
  ::glVertex3d(+1, +1, +1);
  ::glTexCoord2d(0.0, 1.0);
  ::glVertex3d(-1, +1, +1);
  ::glEnd();
  ::glBindTexture(GL_TEXTURE_2D, 0);
  ::glDisable(GL_TEXTURE_2D);
  //
  ::glPopMatrix();
}

void delfem2::opengl::CDrawerOldGL_Render2Tex::Draw(const CRender2Tex &r2t) const {

  ::glPointSize(this->pointSize);
  if (isDrawDepth) {
    this->Draw_Point(r2t);
  }
  // -----------
  ::glLineWidth(3);
  this->Draw_Axis(r2t);
  // ----------
  ::glLineWidth(1);
  ::glColor3d(0, 0, 0);
  delfem2::opengl::CDrawerOldGL_Render2Tex::Draw_BoundingBox(r2t);
  // -----------
  if (r2t.id_tex_color > 0 && this->isDrawTex) {
    ::glBindTexture(GL_TEXTURE_2D, r2t.id_tex_color);
    this->Draw_Texture(r2t);
  }

}

void delfem2::opengl::CDrawerOldGL_Render2Tex::Draw_Axis(const CRender2Tex &r2t) const {
  double mMVP[16];
  MatMat4(
      mMVP,
      r2t.mat_projection,
      r2t.mat_modelview);
  double mMVPinv[16];
  Inverse_Mat4(mMVPinv, mMVP);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glMultMatrixd(TransposeMat4ForOpenGL(mMVPinv,false).data());
  ::glTranslated(-1.01, -1.01, -1.01);
  delfem2::opengl::DrawAxis(draw_len_axis);
  ::glPopMatrix();
}

void delfem2::opengl::CDrawerOldGL_Render2Tex::Draw_BoundingBox(const CRender2Tex &r2t) {
  double mMVP[16];
  MatMat4(
      mMVP,
      r2t.mat_projection,
      r2t.mat_modelview);
  double mMVPinv[16];
  Inverse_Mat4(mMVPinv, mMVP);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glMultMatrixd(TransposeMat4ForOpenGL(mMVPinv,false).data());
  ::glLineWidth(3);
  ::glDisable(GL_LIGHTING);
  const double pmin[3] = {-1., -1., -1.};
  const double pmax[3] = {+1., +1., +1.};
  ::delfem2::opengl::DrawBox3_Edge(pmin, pmax);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPopMatrix();
}

void delfem2::opengl::CDrawerOldGL_Render2Tex::Draw_Point(const CRender2Tex &r2t) const {
  const unsigned int nResX = r2t.width;
  const unsigned int nResY = r2t.height;
  if (r2t.aDepth.size() != nResX * nResY) {
    glPopMatrix();
    return;
  }
  //
  double mMVP[16];
  MatMat4(
      mMVP,
      r2t.mat_projection,
      r2t.mat_modelview);
  double mMVPinv[16];
  Inverse_Mat4(mMVPinv, mMVP);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glMultMatrixd(TransposeMat4ForOpenGL(mMVPinv,false).data());
  ::glDisable(GL_LIGHTING);

  if (colorPoint.size() == 3) { ::glColor3dv(colorPoint.data()); }
  if (colorPoint.size() == 4) { ::glColor4dv(colorPoint.data()); }
  ::glBegin(GL_POINTS);
  for (unsigned int iy = 0; iy < nResY; ++iy) {
    for (unsigned int ix = 0; ix < nResX; ++ix) {
      const double x0 = -1 + 2.0 / nResX * ix;
      const double y0 = -1 + 2.0 / nResY * iy;
      const double z0 = 1 - 2.0 * r2t.aDepth[iy * nResX + ix];
      if (z0 < -0.9 && isDrawOnlyHitPoints) { continue; } // ray is shooted from -1 to +1
      ::glVertex3d(x0, y0, z0);
    }
  }
  ::glEnd();
  // --
  glPopMatrix();
}

// --------------------------------------------------------

void delfem2::opengl::CRender2Tex_DrawOldGL_BOX::Initialize(
    unsigned int nresX_,
    unsigned int nresY_,
    unsigned int nresZ_,
    double elen_) {
  namespace dfm2 = delfem2;
  this->lengrid = elen_;
  //
  aSampler.resize(0);
  aDrawSampler.resize(0);
  {
    aSampler.resize(aSampler.size() + 1);
    auto &smplr = aSampler[aSampler.size() - 1];
    smplr.SetTextureProperty(nresY_, nresZ_, true); // +x
    dfm2::CMat4d::AffineAxisTransform(
        {0,0,1},
        {1,0,0},
        {0,1,0}
    ).CopyTo(smplr.mat_modelview);
    dfm2::CMat4d::AffineOrthogonalProjection(
        -elen_*nresY_*0.5, elen_*nresY_*0.5,
        -elen_*nresZ_*0.5, elen_*nresZ_*0.5,
        -elen_*nresX_*0.5, elen_*nresX_*0.5
    ).CopyTo(smplr.mat_projection);
  }
  {
    aSampler.resize(aSampler.size() + 1);
    auto &smplr = aSampler[aSampler.size() - 1];
    smplr.SetTextureProperty(nresY_, nresZ_, true); // -x
    dfm2::CMat4d::AffineAxisTransform(
        {0,0,-1},
        {1,0,0},
        {0,-1,0}
    ).CopyTo(smplr.mat_modelview);
    dfm2::CMat4d::AffineOrthogonalProjection(
        -elen_*nresY_*0.5, elen_*nresY_*0.5,
        -elen_*nresZ_*0.5, elen_*nresZ_*0.5,
        -elen_*nresX_*0.5, elen_*nresX_*0.5
    ).CopyTo(smplr.mat_projection);
  }
  //
  {
    aSampler.resize(aSampler.size() + 1);
    auto &smplr = aSampler[aSampler.size() - 1];
    smplr.SetTextureProperty(nresX_, nresZ_, true); // +y
    dfm2::CMat4d::AffineAxisTransform(
        {1,0,0},
        {0,0,1},
        {0,-1,0}
    ).CopyTo(smplr.mat_modelview);
    dfm2::CMat4d::AffineOrthogonalProjection(
        -elen_*nresX_*0.5, elen_*nresX_*0.5,
        -elen_*nresZ_*0.5, elen_*nresZ_*0.5,
        -elen_*nresY_*0.5, elen_*nresY_*0.5
    ).CopyTo(smplr.mat_projection);
  }
  {
    aSampler.resize(aSampler.size() + 1);
    auto &smplr = aSampler[aSampler.size() - 1];
    smplr.SetTextureProperty(nresX_, nresZ_, true); // -y
    dfm2::CMat4d::AffineAxisTransform(
        {1,0,0},
        {0,0,-1},
        {0,1,0}
    ).CopyTo(smplr.mat_modelview);
    dfm2::CMat4d::AffineOrthogonalProjection(
        -elen_*nresX_*0.5, elen_*nresX_*0.5,
        -elen_*nresZ_*0.5, elen_*nresZ_*0.5,
        -elen_*nresY_*0.5, elen_*nresY_*0.5
    ).CopyTo(smplr.mat_projection);
  }
  {
    aSampler.resize(aSampler.size() + 1);
    auto &smplr = aSampler[aSampler.size() - 1];
    smplr.SetTextureProperty(nresX_, nresY_, true);
    dfm2::CMat4d::AffineAxisTransform(
        {1,0,0},
        {0,1,0},
        {0,0,1}
    ).CopyTo(smplr.mat_modelview);
    dfm2::CMat4d::AffineOrthogonalProjection(
        -elen_*nresX_*0.5, elen_*nresX_*0.5,
        -elen_*nresY_*0.5, elen_*nresY_*0.5,
        -elen_*nresZ_*0.5, elen_*nresZ_*0.5
    ).CopyTo(smplr.mat_projection);
  }
  {
    aSampler.resize(aSampler.size() + 1);
    auto &smplr = aSampler[aSampler.size() - 1];
    smplr.SetTextureProperty(nresX_, nresY_, true);
    dfm2::CMat4d::AffineAxisTransform(
        {1,0,0},
        {0,-1,0},
        {0,0,-1}
    ).CopyTo(smplr.mat_modelview);
    dfm2::CMat4d::AffineOrthogonalProjection(
        -elen_*nresX_*0.5, elen_*nresX_*0.5,
        -elen_*nresY_*0.5, elen_*nresY_*0.5,
        -elen_*nresZ_*0.5, elen_*nresZ_*0.5
    ).CopyTo(smplr.mat_projection);
  }

  // ------------------------
  aDrawSampler.resize(6);
  aDrawSampler[0].SetPointColor(1.0, 0.0, 0.0);
  aDrawSampler[1].SetPointColor(1.0, 0.5, 0.5);
  aDrawSampler[2].SetPointColor(0.0, 1.0, 0.0);
  aDrawSampler[3].SetPointColor(0.5, 1.0, 0.5);
  aDrawSampler[4].SetPointColor(0.0, 0.0, 1.0);
  aDrawSampler[5].SetPointColor(0.5, 0.5, 1.0);
  for (auto &smplr : aDrawSampler) {
    smplr.draw_len_axis = 0.2;
    smplr.isDrawTex = false;
    smplr.isDrawOnlyHitPoints = true;
  }

}

void delfem2::opengl::CarveVoxelByDepth(
    std::vector<int> &aVal,
    const CRender2Tex_DrawOldGL_BOX &sampler_box) {
  const unsigned int nx = sampler_box.nDivX();
  const unsigned int ny = sampler_box.nDivY();
  const unsigned int nz = sampler_box.nDivZ();
  // ------
  aVal.assign(nz * ny * nx, 1);
  for (unsigned int iy = 0; iy < ny; ++iy) {
    for (unsigned int iz = 0; iz < nz; ++iz) {
      double d0 = sampler_box.aSampler[0].aDepth[iz * ny + iy];
      const unsigned int nd = static_cast<unsigned int>(d0 * nx);
      for (unsigned int id = 0; id < nd; id++) {
        const unsigned int ix0 = nx - id - 1;
        const unsigned int iy0 = iy;
        const unsigned int iz0 = iz;
        aVal[iz0 * ny * nx + iy0 * nx + ix0] = 0;
      }
    }
  }
  for (unsigned int iy = 0; iy < ny; ++iy) {
    for (unsigned int iz = 0; iz < nz; ++iz) {
      double d0 = sampler_box.aSampler[1].aDepth[iz * ny + iy];
      const unsigned int nd = static_cast<unsigned int>(d0 * nx);
      for (unsigned int id = 0; id < nd; id++) {
        const unsigned int ix0 = id;
        const unsigned int iy0 = iy;
        const unsigned int iz0 = nz - 1 - iz;
        aVal[iz0 * ny * nx + iy0 * nx + ix0] = 0;
      }
    }
  }

  for (unsigned int ix = 0; ix < nx; ++ix) {
    for (unsigned int iz = 0; iz < nz; ++iz) {
      double d0 = sampler_box.aSampler[2].aDepth[iz * nx + ix];
      const unsigned int nd = static_cast<unsigned int>(d0 * ny);
      for (unsigned int id = 0; id < nd; id++) {
        const unsigned int ix0 = ix;
        const unsigned int iy0 = ny - 1 - id;
        const unsigned int iz0 = nz - 1 - iz;
        aVal[iz0 * ny * nx + iy0 * nx + ix0] = 0;
      }
    }
  }
  for (unsigned int ix = 0; ix < nx; ++ix) {
    for (unsigned int iz = 0; iz < nz; ++iz) {
      double d0 = sampler_box.aSampler[3].aDepth[iz * nx + ix];
      const unsigned int nd = static_cast<unsigned int>(d0 * ny);
      for (unsigned int id = 0; id < nd; id++) {
        const unsigned int ix0 = ix;
        const unsigned int iy0 = id;
        const unsigned int iz0 = iz;
        aVal[iz0 * ny * nx + iy0 * nx + ix0] = 0;
      }
    }
  }
  for (unsigned int ix = 0; ix < nx; ++ix) {
    for (unsigned int iy = 0; iy < ny; ++iy) {
      double d0 = sampler_box.aSampler[4].aDepth[iy * nx + ix];
      const unsigned int nd = static_cast<unsigned int>(d0 * nz);
      for (unsigned int id = 0; id < nd; id++) {
        const unsigned int ix0 = ix;
        const unsigned int iy0 = iy;
        const unsigned int iz0 = nz - 1 - id;
        aVal[iz0 * ny * nx + iy0 * nx + ix0] = 0;
      }
    }
  }
  for (unsigned int ix = 0; ix < nx; ++ix) {
    for (unsigned int iy = 0; iy < ny; ++iy) {
      double d0 = sampler_box.aSampler[5].aDepth[iy * nx + ix];
      const unsigned int nd = static_cast<unsigned int>(d0 * nz);
      for (unsigned int id = 0; id < nd; id++) {
        const unsigned int ix0 = ix;
        const unsigned int iy0 = ny - 1 - iy;
        const unsigned int iz0 = id;
        aVal[iz0 * ny * nx + iy0 * nx + ix0] = 0;
      }
    }
  }
}
