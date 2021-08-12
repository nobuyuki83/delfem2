#ifndef DFM2_CV_IMGGRAYGRAD_H
#define DFM2_CV_IMGGRAYGRAD_H

namespace delfem2 {

void Grayscale(
    float *gray_,
    int nv,
    int ny,
    int nx,
    const float *aRGB_,
    const float *gray_scale);

void GrayscaleGradientCenter(
    float *grad_,
    int nv,
    int ny,
    int nx,
    const float *aRGB_,
    const float *gray_doeff);

void GrayscaleGradientSobel(
    float *grad_,
    int nv,
    int ny,
    int nx,
    const float *aRGB_,
    const float *gray_coeff);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/cvimggraygrad.cpp"
#endif

#endif