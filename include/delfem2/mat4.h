/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file 3x3 matrix class (CMat3) and functions
 */


#ifndef DFM2_MAT4_H
#define DFM2_MAT4_H

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits> // using NaN Check

#define NEARLY_ZERO 1.e-16

// -----------------------------

namespace delfem2 {


void Mat4Vec3(double vo[3],
              const double M[16],
              const double vi[3]);

// -------------------------------

template <typename T>
void MatVec4(T v[4],
    const T A[16],
    const T x[4]);

template <typename T>
void MatMat4(
    T* C,
    const T* A,
    const T* B);


template <typename REAL>
void Mat4_Identity(
    REAL A[16]);
  
template <typename REAL>
void Mat4_AffineScale(
    REAL A[16],
    REAL s);
  
template <typename REAL>
void Mat4_AffineTranslation(
    REAL A[16],
    REAL dx, REAL dy, REAL dz);

template <typename T>
void Mat4_AffineRotationRodriguez(
    T A[16],
    T dx, T dy, T dz);

/**
 * @func multiply rotation affine matrix from left to an affine matrix in 3D
 * @details the ritation is parmeterized with a rodriguez angle
 */
template <typename REAL>
void Rotate_Mat4AffineRodriguez(
    REAL A[16],
    const REAL v[3]);

/**
 * @func multiply translation affine matrix from left to an affine matrix in 3D
 */
template <typename REAL>
void Translate_Mat4Affine(
    REAL A[16],
    const REAL v[3]);

template <typename T>
void Vec3_Mat4Vec3_AffineProjection(
                                T y0[3],
                                const T a[16],
                                const T x0[3]);

}

#endif
