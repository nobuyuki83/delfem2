/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef MAT3_H
#define MAT3_H

#include <vector>
#include <cassert>
#include <math.h>
#include <iostream>
#include <limits> // using NaN Check

//#include "vec3.h"

#define NEARLY_ZERO 1.e-16

void MatMat3(double* UL,
             const double* U, const double* L);
void MatMatTrans3(double* ULUt,
                  const double* UL, const double* U);

// compute eigen value & vector for symmmetric matrix
// sm[6] = (M_00,M_11,M_22,M_12,M_20,M_01)
// M = ULU^T
// u[9] = (U_00,U_01,U_02, U_10,U_11,U_12, U_20,U_21,U_22)
bool eigenSym3(double u[9], double l[3],
               const double sm[6],
               int nitr);

void svd3(double U[9], double G[3], double V[9],
          const double m[9],
          int nitr);

double Det(const double U[9]);

void GetRotPolarDecomp(double R[9],
                       const double am[9],
                       int nitr);



////////////////////////////////////////////////////////////////////////////

class CMatrix3; // this pre-definition is needed for following functions
CMatrix3 operator+ (const CMatrix3& lhs, const CMatrix3& rhs);
CMatrix3 operator- (const CMatrix3& lhs, const CMatrix3& rhs);
CMatrix3 operator* (double d, const CMatrix3& rhs);
CMatrix3 operator* (const CMatrix3& m, double d);
CMatrix3 operator* (const CMatrix3& lhs, const CMatrix3& rhs);
CMatrix3 operator/ (const CMatrix3& m, double d);
std::ostream &operator<<(std::ostream &output, const CMatrix3& m);
std::istream &operator>>(std::istream &output, CMatrix3& m);

static bool myIsNAN_Matrix3(double d){
  return !(d > d-1);
}
//!< 3x3 matrix
class CMatrix3
{
public:
  CMatrix3();
  CMatrix3(const double s);
  CMatrix3(double v00, double v01, double v02,
           double v10, double v11, double v12,
           double v20, double v21, double v22);
  CMatrix3(double x, double y, double z);
  CMatrix3(const double m[9]);
  ////////
  void GetElements(double m[9]) const { for(unsigned int i=0;i<9;i++){ m[i]=mat[i]; } }
  void AffineMatrixTrans(double m[16]) const {
    m[0*4+0] = mat[0];  m[1*4+0] = mat[1];  m[2*4+0] = mat[2];  m[3*4+0] = 0;
    m[0*4+1] = mat[3];  m[1*4+1] = mat[4];  m[2*4+1] = mat[5];  m[3*4+1] = 0;
    m[0*4+2] = mat[6];  m[1*4+2] = mat[7];  m[2*4+2] = mat[8];  m[3*4+2] = 0;
    m[0*4+3] = 0;       m[1*4+3] = 0;       m[2*4+3] = 0;       m[3*4+3] = 1;
  }
  double Get(int i, int j) const {
    return mat[i*3+j];
  }
  ///////////////////////////////////////////////
//  CVector3 MatVec(const CVector3& vec0) const;
  void MatVec(const double vec0[], double vec1[]) const;
  void MatVecTrans(const double vec0[], double vec1[]) const;
//  CVector3 MatVecTrans(const CVector3& vec0) const;
  CMatrix3 MatMat(const CMatrix3& mat0) const;
  CMatrix3 MatMatTrans(const CMatrix3& mat0) const;
  /////
  CMatrix3 Sym() const{
    CMatrix3 m;
    for(unsigned int i=0;i<3;i++){
    for(unsigned int j=0;j<3;j++){
      m.mat[i*3+j] = (mat[i*3+j]+mat[j*3+i])*0.5;
    }
    }
    return m;
  }
  inline const CMatrix3 operator-() const{ return (*this)*(-1.0); }
  inline const CMatrix3 operator+() const{ return (*this); }
  inline CMatrix3& operator+=(const CMatrix3& rhs){
    for(unsigned int i=0;i<9;i++){ mat[i] += rhs.mat[i]; }
		return *this;
	}  
  inline CMatrix3& operator-=(const CMatrix3& rhs){
    for(unsigned int i=0;i<9;i++){ mat[i] -= rhs.mat[i]; }
		return *this;
	}    
	inline CMatrix3& operator*=(double d){
    for(unsigned int i=0;i<9;i++){ mat[i] *= d; }
		return *this;
	}
	inline CMatrix3& operator/=(double d){
    double invd = 1.0/d;
    for(unsigned int i=0;i<9;i++){ mat[i] *= invd; }
		return *this;
	}
  inline double operator[](int i) const{
    return this->mat[i];
  }  
  inline double& operator()(int i, int j){
    return this->mat[i*3+j];
  }
  /////
  /*
  CVector3 getColumnVector(int j) const {
    return CVector3(mat[0*3+j], mat[1*3+j], mat[2*3+j]);
  }
   */
  CMatrix3 Inverse() const;
  /////
  void SetInverse();
//  void SetDiag(const CVector3& d);
  void SetSymetric(const double sm[6]);
  void SetZero();
  void SetRandom();
//  void SetRotMatrix_Cartesian(const CVector3& v);
  void SetRotMatrix_Cartesian(const double vec[]);
  void SetRotMatrix_Cartesian(double x, double y, double z);
  void SetRotMatrix_Rodrigues(const double vec[]);
  void SetRotMatrix_CRV(const double crv[]);
  void SetRotMatrix_Quaternion(const double quat[]);
  void SetRotMatrix_BryantAngle(double rx, double ry, double rz);
//  void SetSpinTensor(const CVector3& vec0);
//  void SetOuterProduct(const CVector3& vec0, const CVector3& vec1 );
//  void SetProjection(const CVector3& vec0);
  void SetIdentity(double scale = 1);
  /////
//  CVector3 GetCartesianRotationVector() const;
//  CVector3 GetSpinVector()const;
  void GetCRV_RotMatrix(double crv[]) const;
  void GetQuat_RotMatrix(double quat[]) const;
  ///////
  CMatrix3 Trans() const {
    CMatrix3 m;
    m.mat[0] = mat[0]; m.mat[1] = mat[3]; m.mat[2] = mat[6];
    m.mat[3] = mat[1]; m.mat[4] = mat[4]; m.mat[5] = mat[7];
    m.mat[6] = mat[2]; m.mat[7] = mat[5]; m.mat[8] = mat[8];    
    return m;
  }
  bool isNaN() const{
    double s=mat[0]+mat[1]+mat[2]+mat[3]+mat[4]+mat[5]+mat[6]+mat[7]+mat[8];
    return myIsNAN_Matrix3( s ) != 0;
  }
  double Det() const {
    return
    + mat[0]*mat[4]*mat[8] + mat[3]*mat[7]*mat[2] + mat[6]*mat[1]*mat[5]
    - mat[0]*mat[7]*mat[5] - mat[6]*mat[4]*mat[2] - mat[3]*mat[1]*mat[8];
  }
  double SqNorm_Frobenius() const {
    double s = 0.0;
    for(int i=0;i<9;++i){
      s += mat[i]*mat[i];
    }
    return s;
  }
  double Trace() const {
    return mat[0]+mat[4]+mat[8];
  }
  double SecondInvarint() const {
    const CMatrix3& m2 = (*this)*(*this);
    const double tr = this->Trace();
    return 0.5*(tr*tr-m2.Trace());
  }
  void Print() const {
    std::cout << mat[0] << " " << mat[1] << " " << mat[2] << std::endl;
    std::cout << mat[3] << " " << mat[4] << " " << mat[5] << std::endl;
    std::cout << mat[6] << " " << mat[7] << " " << mat[8] << std::endl;
  }
  //////////////////////////////////////
  // static functions
  static CMatrix3 Identity(double scale = 1){
    CMatrix3 m;
    m.SetIdentity(scale);
    return m;
  }
  static CMatrix3 Zero(){
    CMatrix3 m;
    m.SetZero();
    return m;
  }

public:
  double mat[9];  
};


/////////////////////////////////////////////////////////////////////

void MatVec4(double v[4],
             const double A[16],
             const double x[4]);
void Affine3D(double y0[3],
              const double a[16],
              const double x0[3]);
void SetAffine_Scale(double A[16],
                     double s);
void SetAffine_Trans(double A[16],
                     double dx, double dy, double dz);
void SetAffine_Rotate_Rodriguez(double A[16],
                                double dx, double dy, double dz);

#endif
