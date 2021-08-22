/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// Done: having MVP matrix instead of axes

#ifndef DFM2_OPENGL_R2T_H
#define DFM2_OPENGL_R2T_H

#include <cstdio>
#include <vector>

#include "delfem2/mat4.h"
#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2 {

DFM2_INLINE void Mat4_OrthongoalProjection_AffineTrans(
    double mMV[16],
    double mP[16],
    const double origin[3],
    const double az[3],
    const double ax[3],
    unsigned int nResX,
    unsigned int nResY,
    double lengrid,
    double z_range);

namespace opengl {

/**
 * Class to prepare framebuffer object (FBO) for "render-to-texture" purpose
 * The rendered image and depth are saved in the CPU memory at the end.
 */
class CRender2Tex {
 public:
  CRender2Tex() :
      mat_modelview_colmajor{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
      mat_projection_colmajor{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1} {
    width = 128;
    height = 128;
    is_rgba_8ui = true;
    id_tex_color = 0;
    id_tex_depth = 0;
    id_framebuffer = 0;
  }
  // --------------------------
  void InitGL();
  void SetZeroToDepth() { for (float &i : aZ) { i = 0.0; }}
  void GetMVPG(double mMVPG[16]) const {
    double mMVP[16];
    MatMat4(mMVP, mat_modelview_colmajor, mat_projection_colmajor);
    const double tmp0 = width * 0.5;
    const double tmp1 = height * 0.5;
    double mG[16] = {
        tmp0, 0, 0, 0,
        0, tmp1, 0, 0,
        0, 0, 0.5, 0,
        tmp0 - 0.5, tmp1 - 0.5, 0.5, 1};
    MatMat4(mMVPG, mMVP, mG);
  }
  /**
  * @brief update the bounding box by adding points
  * @param[in,out] pmin lower coner
  * @param[in,out] pmax upper corner
  * @details if( pmin[0] > pmax[0] ) this bounding box is empty
  */
  DFM2_INLINE void BoundingBox3(double *pmin, double *pmax) const;

  void SetTextureProperty(unsigned int nw, unsigned int nh, bool is_rgba_8ui_) {
    this->width = nw;
    this->height = nh;
    this->is_rgba_8ui = is_rgba_8ui_;
  }
  void SetValue_CpuImage_8ui(
      const unsigned char *image,
      unsigned int nw,
      unsigned int nh,
      unsigned int nch) {
    aRGBA_8ui.resize(nw * nh * 4, 255);
    for (unsigned int ih = 0; ih < nh; ++ih) {
      for (unsigned int iw = 0; iw < nw; ++iw) {
        aRGBA_8ui[(ih * nw + iw) * 4 + 0] = image[(ih * nw + iw) * 3 + 0];
        aRGBA_8ui[(ih * nw + iw) * 4 + 1] = image[(ih * nw + iw) * 3 + 1];
        aRGBA_8ui[(ih * nw + iw) * 4 + 2] = image[(ih * nw + iw) * 3 + 2];
      }
    }
  }
  void Start();
  void End();
  template<typename REAL>
  std::vector<REAL> GetAffineMatrixModelViewAsColMajorStlVector() const {
    return std::vector<REAL>(mat_modelview_colmajor, mat_modelview_colmajor + 16);
  }
  template<typename REAL>
  std::vector<REAL> GetAffineMatrixProjectionAsColMajorStlVector() const {
    return std::vector<REAL>(mat_projection_colmajor, mat_projection_colmajor + 16);
  }
  template<typename REAL>
  void SetAffineMatrixModelViewAsColMajorStlVector(const std::vector<REAL>& a) {
    assert(a.size()==16);
    for(int i=0;i<16;++i){ mat_modelview_colmajor[i] = a[i]; }
  }
  template<typename REAL>
  void SetAffineMatrixProjectionAsColMajorStlVector(const std::vector<REAL>& a) {
    assert(a.size()==16);
    for(int i=0;i<16;++i){ mat_projection_colmajor[i] = a[i]; }
  }
 private:
  void CopyToCPU_Depth();
  void CopyToCPU_RGBA8UI();
  void CopyToCPU_RGBA32F();
 public:
  unsigned int width;
  unsigned int height;
  // ------------------
  unsigned int id_tex_color;
  unsigned int id_tex_depth;
  unsigned int id_framebuffer;
  // --------
  std::vector<float> aZ;
  bool is_rgba_8ui;
  std::vector<unsigned char> aRGBA_8ui;
  std::vector<float> aRGBA_32f;
  //
  double mat_modelview_colmajor[16];  // affine matrix
  double mat_projection_colmajor[16];  // affine matrix
//protected:
  int view[4]{}; // viewport information
};

/**
 * @brief project input point to the depth surface
 * @param[in] ps the point to project
 */
DFM2_INLINE bool GetProjectedPoint(
    CVec3d &p0,
    CVec3d &n0,
    const CVec3d &ps,
    const CRender2Tex &smplr);

} // opengl
} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/r2t.cpp"
#endif

#endif
