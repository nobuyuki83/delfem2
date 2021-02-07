
#include "delfem2/cvimggraygrad.h"
#include "delfem2/viewtensor.h"
#include "delfem2/dfm2_inline.h"
#include <cassert>
#include <cmath>

namespace delfem2 {
namespace cvimggraygrad {

inline float Dot3(const float p[3], const float q[3]) {
  return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
}

}
}

DFM2_INLINE void delfem2::Grayscale(
    float* gray_,
    int nv,
    int ny,
    int nx,
    const float* aRGB_,
    const float *gray_coeff)
{
  delfem2::ViewTensor3<float> aGray(gray_, nv, ny, nx);
  const delfem2::ViewTensor4Const<float> frame(aRGB_, nv, ny, nx, 3);

  for (int iv = 0; iv < nv; ++iv) {
    for (int iy1 = 0; iy1 < ny; ++iy1) {
      for (int ix1 = 0; ix1 < nx; ++ix1) {
        aGray(iv,iy1,ix1) = cvimggraygrad::Dot3(frame.data(iv, iy1, ix1), gray_coeff);
      }
    }
  }
}


DFM2_INLINE void delfem2::GrayscaleGradientCenter(
    float* grad_,
    int nv,
    int ny,
    int nx,
    const float* aRGB_,
    const float *gray_coeff)
{
  delfem2::ViewTensor4<float> aGrad(grad_,nv,ny,nx,2);
  const delfem2::ViewTensor4Const<float> frame(aRGB_, nv, ny, nx, 3);

  for(int iv = 0; iv < nv; ++iv){
    for (int iy1 = 0; iy1 < ny; ++iy1) {
      for (int ix1 = 0; ix1 < nx; ++ix1) {
        const int ix0 = (ix1 == 0) ? 0 : ix1-1; // reflection
        const int iy0 = (iy1 == 0) ? 0 : iy1-1; // reflection
        const int ix2 = (ix1 == nx - 1) ? nx - 1 : ix1 + 1;
        const int iy2 = (iy1 == ny - 1) ? ny - 1 : iy1 + 1;
        assert( ix0 >=0 && ix0 < nx );
        assert( iy0 >=0 && iy0 < ny );
        assert( ix1 >=0 && ix1 < nx );
        assert( iy1 >=0 && iy1 < ny );
        const auto v01 = cvimggraygrad::Dot3(frame.data(iv,iy0,ix1),gray_coeff);
        const auto v10 = cvimggraygrad::Dot3(frame.data(iv,iy1,ix0),gray_coeff);
        const auto v12 = cvimggraygrad::Dot3(frame.data(iv,iy1,ix2),gray_coeff);
        const auto v21 = cvimggraygrad::Dot3(frame.data(iv,iy2,ix1),gray_coeff);
        aGrad(iv, iy1, ix1, 0) = (v12 - v10)*0.5;
        aGrad(iv, iy1, ix1, 1) = (v21 - v01)*0.5;
      }
    }
  }
}

DFM2_INLINE void delfem2::GrayscaleGradientSobel(
    float* grad_,
    int nv,
    int ny,
    int nx,
    const float* aRGB_,
    const float gray_coeff[3])
{
  delfem2::ViewTensor4<float> aGrad(grad_,nv,ny,nx,2);
  const delfem2::ViewTensor4Const<float> frame(aRGB_, nv, ny, nx, 3);

  for(int iv = 0; iv < nv; ++iv){
    for (int iy1 = 0; iy1 < ny; ++iy1) {
      for (int ix1 = 0; ix1 < nx; ++ix1) {
        const int ix0 = abs(ix1-1); // reflection
        const int iy0 = abs(iy1-1); // reflection
        const int ix2 = nx - 1 - ::abs(nx - 2 - ix1); // reflection
        const int iy2 = ny - 1 - ::abs(ny - 2 - iy1); // reflection
        assert( ix0 >=0 && ix0 < nx );
        assert( iy0 >=0 && iy0 < ny );
        assert( ix1 >=0 && ix1 < nx );
        assert( iy1 >=0 && iy1 < ny );
        const auto v00 = cvimggraygrad::Dot3(frame.data(iv,iy0,ix0),gray_coeff);
        const auto v01 = cvimggraygrad::Dot3(frame.data(iv,iy0,ix1),gray_coeff);
        const auto v02 = cvimggraygrad::Dot3(frame.data(iv,iy0,ix2),gray_coeff);
        const auto v10 = cvimggraygrad::Dot3(frame.data(iv,iy1,ix0),gray_coeff);
        const auto v12 = cvimggraygrad::Dot3(frame.data(iv,iy1,ix2),gray_coeff);
        const auto v20 = cvimggraygrad::Dot3(frame.data(iv,iy2,ix0),gray_coeff);
        const auto v21 = cvimggraygrad::Dot3(frame.data(iv,iy2,ix1),gray_coeff);
        const auto v22 = cvimggraygrad::Dot3(frame.data(iv,iy2,ix2),gray_coeff);
        const float dx = (-v00 - 2*v10 - v20 + v02 + 2*v12 + v22)/8.f;
        const float dy = (-v00 - 2*v01 - v02 + v20 + 2*v21 + v22)/8.f;
        aGrad(iv, iy1, ix1, 0) = dx;
        aGrad(iv, iy1, ix1, 1) = dy;
      }
    }
  }
}