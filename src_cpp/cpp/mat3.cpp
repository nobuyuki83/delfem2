/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>

#include "delfem2/mat3.h"

void MatMat3
(double* C,
 const double* A, const double* B)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i*3+j] = A[i*3+0]*B[0*3+j] + A[i*3+1]*B[1*3+j] + A[i*3+2]*B[2*3+j];
    }
  }
}

void MatMatTrans3
(double* C,
 const double* A, const double* B)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i*3+j] = A[i*3+0]*B[j*3+0] + A[i*3+1]*B[j*3+1] + A[i*3+2]*B[j*3+2];
    }
  }
}

double sqnormFrobenius
(const double sm[6])
{
  return sm[0]*sm[0] + sm[1]*sm[1] + sm[2]*sm[2] + 2*(sm[3]*sm[3] + sm[4]*sm[4] + sm[5]*sm[5]);
}

// compute eigen value & vector for symmmetric matrix
// sm[6] = (M_00,M_11,M_22,M_12,M_20,M_01)
// M = ULU^T
// u[9] = (U_00,U_01,U_02, U_10,U_11,U_12, U_20,U_21,U_22)
bool eigenSym3
(double u[9], double l[3],
 const double sm[6],
 int nitr)
{
  //  double tol = 1.0e-5;
  u[0]=u[4]=u[8]=1.0;
  u[1]=u[2]=u[3]=u[5]=u[6]=u[7]=0.0;
  l[0]=l[1]=l[2]=0.0;
  double dnrm = sqnormFrobenius(sm);
  if( dnrm < 1.0e-20 ){
    // this matrix is too small;
    return false;
  }
  double scale = sqrt(dnrm);
  double invscl = 1.0/scale;
  double sms[6] = {sm[0]*invscl,sm[1]*invscl,sm[2]*invscl,sm[3]*invscl,sm[4]*invscl,sm[5]*invscl};
  for(int itr=0;itr<nitr;itr++){
    const double m[6] = {sms[0],sms[1],sms[2],sms[3],sms[4],sms[5]};
    const double v[9] = {u[0],u[1],u[2], u[3],u[4],u[5], u[6],u[7],u[8] };
    const double a12 = fabs(sms[3]);
    const double a20 = fabs(sms[4]);
    const double a01 = fabs(sms[5]);
    if( a12 > a20 && a12 > a01 ){ // a12 sms[3] is biggest
      const double t = 0.5*atan2(2*m[3],m[2]-m[1]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[1] = ct*ct*m[1]+st*st*m[2]-2*st*ct*m[3];
      sms[2] = ct*ct*m[2]+st*st*m[1]+2*st*ct*m[3];
      sms[3] = 0; // (ct*ct-st*st)*m[3]+st*ct*(m[1]-m[2]);
      sms[4] = st*m[5]+ct*m[4];
      sms[5] = ct*m[5]-st*m[4];
      ////
      u[1] = +ct*v[1]-st*v[2];
      u[2] = +st*v[1]+ct*v[2];
      u[4] = +ct*v[4]-st*v[5];
      u[5] = +st*v[4]+ct*v[5];
      u[7] = +ct*v[7]-st*v[8];
      u[8] = +st*v[7]+ct*v[8];
    }
    if( a20 > a01 && a20 > a12 ){ // a20 sms[4] is biggest
      const double t = 0.5*atan2(2*m[4],m[2]-m[0]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[0] = ct*ct*m[0]+st*st*m[2]-2*st*ct*m[4];
      sms[2] = ct*ct*m[2]+st*st*m[0]+2*st*ct*m[4];
      sms[3] = st*m[5]+ct*m[3];
      sms[4] = 0; // (ct*ct-st*st)*m[4]+st*ct*(m[0]-m[2]);
      sms[5] = ct*m[5]-st*m[3];
      ////
      u[0] = +ct*v[0]-st*v[2];
      u[2] = +st*v[0]+ct*v[2];
      u[3] = +ct*v[3]-st*v[5];
      u[5] = +st*v[3]+ct*v[5];
      u[6] = +ct*v[6]-st*v[8];
      u[8] = +st*v[6]+ct*v[8];
    }
    if( a01 > a12 && a01 > a20 ){ // a01 sms[5] is biggest
      const double t = 0.5*atan2(2*m[5],m[1]-m[0]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[0] = ct*ct*m[0]+st*st*m[1]-2*st*ct*m[5];
      sms[1] = ct*ct*m[1]+st*st*m[0]+2*st*ct*m[5];
      sms[3] = st*m[4]+ct*m[3];
      sms[4] = ct*m[4]-st*m[3];
      sms[5] = 0; // (ct*ct-st*st)*m[5]+st*ct*(m[0]-m[1]);
      /////
      u[0] = +ct*v[0]-st*v[1];
      u[1] = +st*v[0]+ct*v[1];
      u[3] = +ct*v[3]-st*v[4];
      u[4] = +st*v[3]+ct*v[4];
      u[6] = +ct*v[6]-st*v[7];
      u[7] = +st*v[6]+ct*v[7];
    }
  }
  l[0] = scale*sms[0];
  l[1] = scale*sms[1];
  l[2] = scale*sms[2];
  return true;
}

double estimationMaxEigenValue(const double mtm[6])
{
  double maxl = 1;
  {  // estimation of maximum eigen value using Gerschgorin's circle theorem
    maxl = mtm[0] + fabs(mtm[3])+fabs(mtm[5]);
    const double tmp2 = mtm[1] + fabs(mtm[3])+fabs(mtm[4]);
    maxl = ( tmp2 > maxl ) ? tmp2 : maxl;
    const double tmp3 = mtm[2] + fabs(mtm[5])+fabs(mtm[4]);
    maxl = ( tmp3 > maxl ) ? tmp3 : maxl;
  }
  return maxl;
}

double Det(const double U[9])
{
  return + U[0]*U[4]*U[8] + U[3]*U[7]*U[2] + U[6]*U[1]*U[5]
  - U[0]*U[7]*U[5] - U[6]*U[4]*U[2] - U[3]*U[1]*U[8];
}

void svd3
(double U[9], double G[3], double V[9],
 const double m[9],
 int nitr)
{
  const double mtm[6] = {
    m[0]*m[0]+m[3]*m[3]+m[6]*m[6],
    m[1]*m[1]+m[4]*m[4]+m[7]*m[7],
    m[2]*m[2]+m[5]*m[5]+m[8]*m[8],
    m[1]*m[2]+m[4]*m[5]+m[7]*m[8],
    m[2]*m[0]+m[5]*m[3]+m[8]*m[6],
    m[0]*m[1]+m[3]*m[4]+m[6]*m[7] };
  double lv[3]; eigenSym3(V,lv,
                          mtm,nitr);
  G[0] = sqrt(lv[0]);
  G[1] = sqrt(lv[1]);
  G[2] = sqrt(lv[2]);
  if( G[0]>1.e-20 && G[1]>1.e-20 && G[2]>1.e-20 ){
    const double invG[3] = { 1.0/G[0], 1.0/G[1], 1.0/G[2] };
    U[0] = (m[0]*V[0]+m[1]*V[3]+m[2]*V[6])*invG[0];
    U[1] = (m[0]*V[1]+m[1]*V[4]+m[2]*V[7])*invG[1];
    U[2] = (m[0]*V[2]+m[1]*V[5]+m[2]*V[8])*invG[2];
    U[3] = (m[3]*V[0]+m[4]*V[3]+m[5]*V[6])*invG[0];
    U[4] = (m[3]*V[1]+m[4]*V[4]+m[5]*V[7])*invG[1];
    U[5] = (m[3]*V[2]+m[4]*V[5]+m[5]*V[8])*invG[2];
    U[6] = (m[6]*V[0]+m[7]*V[3]+m[8]*V[6])*invG[0];
    U[7] = (m[6]*V[1]+m[7]*V[4]+m[8]*V[7])*invG[1];
    U[8] = (m[6]*V[2]+m[7]*V[5]+m[8]*V[8])*invG[2];
  }
  else{
    const double mmt[6] = {
      m[0]*m[0]+m[3]*m[3]+m[6]*m[6],
      m[1]*m[1]+m[4]*m[4]+m[7]*m[7],
      m[2]*m[2]+m[5]*m[5]+m[8]*m[8],
      m[3]*m[6]+m[4]*m[7]+m[5]*m[8],
      m[6]*m[0]+m[7]*m[1]+m[8]*m[2],
      m[0]*m[3]+m[1]*m[4]+m[2]*m[5] };
    double lu[3]; eigenSym3(U,lu,
                            mmt,nitr);
  }
  const double detU = Det(U);
  if( detU < 0 ){
    U[2]=-U[2];
    U[5]=-U[5];
    U[8]=-U[8];
    G[2]=-G[2];
  }
}

void GetRotPolarDecomp
(double R[9],
 ////
 const double am[9],
 int nitr)
{
  double U[9], G[3], V[9];
  svd3(U,G,V,
       am,nitr);
  MatMatTrans3(R,U,V);
}

void SetMatrix3_Quaternion(double r[], const double q[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  r[0] = 1.0 - y2 - z2;
  r[1] = xy - zw;
  r[2] = zx + yw;
  r[3] = xy + zw;
  r[4] = 1.0 - z2 - x2;
  r[5] = yz - xw;
  r[6] = zx - yw;
  r[7] = yz + xw;
  r[8] = 1.0 - x2 - y2;
}

//////////////////////////////////////////////////////////////

CMatrix3 operator* (double d, const CMatrix3& rhs){
  CMatrix3 temp = rhs;
  temp *= d;
  return temp;
}

CMatrix3 operator* (const CMatrix3& m, double d){
  CMatrix3 t = m;
  t *= d;
  return t;
}

CMatrix3 operator/ (const CMatrix3& m, double d){
  CMatrix3 temp = m;
  temp /= d;
  return temp;
}

CMatrix3 operator+ (const CMatrix3& lhs, const CMatrix3& rhs){
  CMatrix3 temp = lhs;
  temp += rhs;
  return temp;
}

CMatrix3 operator* (const CMatrix3& lhs, const CMatrix3& rhs){
  return lhs.MatMat(rhs);
}

CMatrix3 operator- (const CMatrix3& lhs, const CMatrix3& rhs){
  CMatrix3 temp = lhs;
  temp -= rhs;
  return temp;
}

std::ostream &operator<<(std::ostream &output, const CMatrix3& m)
{
  output.setf(std::ios::scientific);
  output << m.mat[0*3+0] << " " << m.mat[0*3+1] << " " << m.mat[0*3+2] << " ";
  output << m.mat[1*3+0] << " " << m.mat[1*3+1] << " " << m.mat[1*3+2] << " ";
  output << m.mat[2*3+0] << " " << m.mat[2*3+1] << " " << m.mat[2*3+2] << " ";
  return output;
}

std::istream &operator>>(std::istream &input, CMatrix3& m)
{
  input >> m.mat[0*3+0] >> m.mat[0*3+1] >> m.mat[0*3+2];
  input >> m.mat[1*3+0] >> m.mat[1*3+1] >> m.mat[1*3+2];
  input >> m.mat[2*3+0] >> m.mat[2*3+1] >> m.mat[2*3+2];
  return input;
}

//////////

CMatrix3::CMatrix3(){
  mat[0*3+0]=0; mat[0*3+1]=0; mat[0*3+2]=0;
  mat[1*3+0]=0; mat[1*3+1]=0; mat[1*3+2]=0;
  mat[2*3+0]=0; mat[2*3+1]=0; mat[2*3+2]=0;
}
CMatrix3::CMatrix3(const double s){
  mat[0*3+0]=s; mat[0*3+1]=0; mat[0*3+2]=0;
  mat[1*3+0]=0; mat[1*3+1]=s; mat[1*3+2]=0;
  mat[2*3+0]=0; mat[2*3+1]=0; mat[2*3+2]=s;
}
CMatrix3::CMatrix3(double v00, double v01, double v02,
                   double v10, double v11, double v12,
                   double v20, double v21, double v22)
{
  mat[0*3+0]=v00; mat[0*3+1]=v01; mat[0*3+2]=v02;
  mat[1*3+0]=v10; mat[1*3+1]=v11; mat[1*3+2]=v12;
  mat[2*3+0]=v20; mat[2*3+1]=v21; mat[2*3+2]=v22;
}

CMatrix3::CMatrix3(double x, double y, double z){
  mat[0*3+0]=x; mat[0*3+1]=0; mat[0*3+2]=0;
  mat[1*3+0]=0; mat[1*3+1]=y; mat[1*3+2]=0;
  mat[2*3+0]=0; mat[2*3+1]=0; mat[2*3+2]=z;
}
CMatrix3::CMatrix3(const double m[9]){
  for(unsigned int i=0;i<9;i++){ mat[i]=m[i]; }
}





void CMatrix3::MatVec(const double vec0[], double vec1[]) const
{
  vec1[0] = mat[0]*vec0[0] + mat[1]*vec0[1] + mat[2]*vec0[2];
  vec1[1] = mat[3]*vec0[0] + mat[4]*vec0[1] + mat[5]*vec0[2];
  vec1[2] = mat[6]*vec0[0] + mat[7]*vec0[1] + mat[8]*vec0[2];
}

void CMatrix3::MatVecTrans(const double vec0[], double vec1[]) const
{
  vec1[0] = mat[0]*vec0[0] + mat[3]*vec0[1] + mat[6]*vec0[2];
  vec1[1] = mat[1]*vec0[0] + mat[4]*vec0[1] + mat[7]*vec0[2];
  vec1[2] = mat[2]*vec0[0] + mat[5]*vec0[1] + mat[8]*vec0[2];
}





CMatrix3 CMatrix3::MatMat(const CMatrix3& mat0) const{
  CMatrix3 m;
  for(unsigned int i=0;i<3;i++){
    for(unsigned int j=0;j<3;j++){
      m.mat[i*3+j] =
      mat[i*3+0]*mat0.mat[0*3+j]
      + mat[i*3+1]*mat0.mat[1*3+j]
      + mat[i*3+2]*mat0.mat[2*3+j];
    }
  }
  return m;
}

CMatrix3 CMatrix3::MatMatTrans(const CMatrix3& mat0) const
{
  CMatrix3 m;
  for(unsigned int i=0;i<3;i++){
    for(unsigned int j=0;j<3;j++){
      m.mat[i*3+j] =
      + mat[0*3+i]*mat0.mat[0*3+j]
      + mat[1*3+i]*mat0.mat[1*3+j]
      + mat[2*3+i]*mat0.mat[2*3+j];
    }
  }
  return m;
}


CMatrix3 CMatrix3::Inverse() const
{
  CMatrix3 mi = *this;
  mi.SetInverse();
  return mi;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void CMatrix3::SetInverse()
{
  const double det = this->Det();
  const double inv_det = 1.0/det;
  double t[9]; for(int i=0;i<9;i++){ t[i] = this->mat[i]; }
  mat[0] = inv_det*(t[4]*t[8]-t[5]*t[7]);
  mat[1] = inv_det*(t[2]*t[7]-t[1]*t[8]);
  mat[2] = inv_det*(t[1]*t[5]-t[2]*t[4]);
  
  mat[3] = inv_det*(t[5]*t[6]-t[3]*t[8]);
  mat[4] = inv_det*(t[0]*t[8]-t[2]*t[6]);
  mat[5] = inv_det*(t[2]*t[3]-t[0]*t[5]);
  
  mat[6] = inv_det*(t[3]*t[7]-t[4]*t[6]);
  mat[7] = inv_det*(t[1]*t[6]-t[0]*t[7]);
  mat[8] = inv_det*(t[0]*t[4]-t[1]*t[3]);
}



void CMatrix3::SetSymetric(const double sm[6])
{
  mat[0] = sm[0];
  mat[1] = sm[5];
  mat[2] = sm[4];
  mat[3] = sm[5];
  mat[4] = sm[1];
  mat[5] = sm[3];
  mat[6] = sm[4];
  mat[7] = sm[3];
  mat[8] = sm[2];
}

void CMatrix3::SetZero()
{
  for(int i=0;i<9;++i){ mat[i] = 0.0; }
}


void CMatrix3::SetRandom(){
  for(int i=0;i<9;i++){
    mat[i] = ((double)std::rand()/(RAND_MAX+1.0))*100-50;
  }
}

void CMatrix3::SetRotMatrix_Cartesian(const double vec[])
{
  double sqt = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
  if( sqt < 1.0e-20 ){ // infinitesmal rotation approximation
    mat[0] = 1;        mat[1] = -vec[2];  mat[2] = +vec[1];
    mat[3] = +vec[2];  mat[4] = 1;        mat[5] = -vec[0];
    mat[6] = -vec[1];  mat[7] = +vec[0];  mat[8] = 1;
    return;
  }
  double t = sqrt(sqt);
  double invt = 1.0/t;
  double n[3] = { vec[0]*invt, vec[1]*invt, vec[2]*invt };
  const double c0 = cos(t);
  const double s0 = sin(t);
  mat[0*3+0] = c0        +(1-c0)*n[0]*n[0];
  mat[0*3+1] =   -n[2]*s0+(1-c0)*n[0]*n[1];
  mat[0*3+2] =   +n[1]*s0+(1-c0)*n[0]*n[2];
  mat[1*3+0] =   +n[2]*s0+(1-c0)*n[1]*n[0];
  mat[1*3+1] = c0        +(1-c0)*n[1]*n[1];
  mat[1*3+2] =   -n[0]*s0+(1-c0)*n[1]*n[2];
  mat[2*3+0] =   -n[1]*s0+(1-c0)*n[2]*n[0];
  mat[2*3+1] =   +n[0]*s0+(1-c0)*n[2]*n[1];
  mat[2*3+2] = c0        +(1-c0)*n[2]*n[2];
}

void CMatrix3::SetRotMatrix_Cartesian(double x, double y, double z){
  const double vec[3] = { x, y, z };
  this->SetRotMatrix_Cartesian(vec);
}



void CMatrix3::SetRotMatrix_Rodrigues(const double vec[])
{
  const double sqlen = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
  const double tmp1 = 1.0/(1+0.25*sqlen);
  mat[0] = 1+tmp1*(       +0.5*vec[0]*vec[0]-0.5*sqlen);
  mat[1] =  +tmp1*(-vec[2]+0.5*vec[0]*vec[1]          );
  mat[2] =  +tmp1*(+vec[1]+0.5*vec[0]*vec[2]          );
  mat[3] =  +tmp1*(+vec[2]+0.5*vec[1]*vec[0]          );
  mat[4] = 1+tmp1*(       +0.5*vec[1]*vec[1]-0.5*sqlen);
  mat[5] =  +tmp1*(-vec[0]+0.5*vec[1]*vec[2]          );
  mat[6] =  +tmp1*(-vec[1]+0.5*vec[2]*vec[0]          );
  mat[7] =  +tmp1*(+vec[0]+0.5*vec[2]*vec[1]          );
  mat[8] = 1+tmp1*(       +0.5*vec[2]*vec[2]-0.5*sqlen);
}

void CMatrix3::SetRotMatrix_CRV(const double crv[])
{
  const double c0 = 0.125*( 16.0 - crv[0]*crv[0] - crv[1]*crv[1] - crv[2]*crv[2] );
  const double tmp = 1.0/( (4.0-c0)*(4.0-c0) );
  mat[0*3+0] = tmp*( (c0*c0+8*c0-16) + 2*crv[0]*crv[0] );
  mat[0*3+1] = tmp*(                   2*crv[0]*crv[1] - 2*c0*crv[2] );
  mat[0*3+2] = tmp*(                   2*crv[0]*crv[2] + 2*c0*crv[1] );
  mat[1*3+0] = tmp*(                   2*crv[1]*crv[0] + 2*c0*crv[2] );
  mat[1*3+1] = tmp*( (c0*c0+8*c0-16) + 2*crv[1]*crv[1] );
  mat[1*3+2] = tmp*(                   2*crv[1]*crv[2] - 2*c0*crv[0] );
  mat[2*3+0] = tmp*(                   2*crv[2]*crv[0] - 2*c0*crv[1] );
  mat[2*3+1] = tmp*(                   2*crv[2]*crv[1] + 2*c0*crv[0] );
  mat[2*3+2] = tmp*( (c0*c0+8*c0-16) + 2*crv[2]*crv[2] );
}

void CMatrix3::SetRotMatrix_Quaternion(const double quat[]){
  SetMatrix3_Quaternion(mat, quat);
}

void CMatrix3::SetRotMatrix_BryantAngle(double rx, double ry, double rz)
{
  CMatrix3 mx; double rvx[3] = {rx,0,0}; mx.SetRotMatrix_Cartesian(rvx);
  CMatrix3 my; double rvy[3] = {0,ry,0}; my.SetRotMatrix_Cartesian(rvy);
  CMatrix3 mz; double rvz[3] = {0,0,rz}; mz.SetRotMatrix_Cartesian(rvz);
  CMatrix3 m = mz;
  m = m.MatMat(my);
  m = m.MatMat(mx);
  *this = m;
}

void CMatrix3::GetQuat_RotMatrix(double quat[]) const{
  const double smat[16] = {
    1+mat[0*3+0]+mat[1*3+1]+mat[2*3+2],
    mat[2*3+1]-mat[1*3+2],
    mat[0*3+2]-mat[2*3+0],
    mat[1*3+0]-mat[0*3+1],
    mat[2*3+1]-mat[1*3+2],
    1+mat[0*3+0]-mat[1*3+1]-mat[2*3+2],
    mat[0*3+1]+mat[1*3+0],
    mat[0*3+2]+mat[2*3+0],
    mat[0*3+2]-mat[2*3+0],
    mat[1*3+0]+mat[0*3+1],
    1-mat[0*3+0]+mat[1*3+1]-mat[2*3+2],
    mat[1*3+2]+mat[2*3+1],
    mat[1*3+0]-mat[0*3+1],
    mat[0*3+2]+mat[2*3+0],
    mat[1*3+2]+mat[2*3+1],
    1-mat[0*3+0]-mat[1*3+1]+mat[2*3+2],
  };
  
  unsigned int imax;
  imax = ( smat[0   *4+   0] > smat[1*4+1] ) ? 0    : 1;
  imax = ( smat[imax*4+imax] > smat[2*4+2] ) ? imax : 2;
  imax = ( smat[imax*4+imax] > smat[3*4+3] ) ? imax : 3;
  
  quat[imax] = 0.5*sqrt(smat[imax*4+imax]);
  for(unsigned int k=0;k<4;k++){
    if( k==imax ) continue;
    quat[k] = smat[imax*4+k]*0.25/quat[imax];
  }
}

void CMatrix3::GetCRV_RotMatrix(double crv[]) const{
  double eparam2[4];
  this->GetQuat_RotMatrix(eparam2);
  crv[0] = 4*eparam2[1]/(1+eparam2[0]);
  crv[1] = 4*eparam2[2]/(1+eparam2[0]);
  crv[2] = 4*eparam2[3]/(1+eparam2[0]);
}


void CMatrix3::SetIdentity(double scale)
{
  mat[0] = scale; mat[1] = 0;     mat[2] = 0;
  mat[3] = 0;     mat[4] = scale; mat[5] = 0;
  mat[6] = 0;     mat[7] = 0;     mat[8] = scale;
}





/////////////////////////////////////////////////////
// 4x4 matrix

void MatVec4
(double v[4],
 const double A[16],
 const double x[4])
{
  v[0] = A[0*4+0]*x[0] + A[0*4+1]*x[1] + A[0*4+2]*x[2] + A[0*4+3]*x[3];
  v[1] = A[1*4+0]*x[0] + A[1*4+1]*x[1] + A[1*4+2]*x[2] + A[1*4+3]*x[3];
  v[2] = A[2*4+0]*x[0] + A[2*4+1]*x[1] + A[2*4+2]*x[2] + A[2*4+3]*x[3];
  v[3] = A[3*4+0]*x[0] + A[3*4+1]*x[1] + A[3*4+2]*x[2] + A[3*4+3]*x[3];
}

void Affine3D
(double y0[3],
 const double a[16],
 const double x0[3])
{
  const double x1[4] = {x0[0], x0[1], x0[2], 1.0};
  double y1[4]; MatVec4(y1,a,x1);
  y0[0] = y1[0]/y1[3];
  y0[1] = y1[1]/y1[3];
  y0[2] = y1[2]/y1[3];
}

void SetAffine_Scale
(double A[16],
 double s)
{
  for(int i=0;i<16;++i){ A[i] = 0.0; }
  A[0*4+0] = s;
  A[1*4+1] = s;
  A[2*4+2] = s;
  A[3*4+3] = 1.0;
}

void SetAffine_Trans
(double A[16],
 double dx, double dy, double dz)
{
  for(int i=0;i<16;++i){ A[i] = 0.0; }
  for(int i=0;i<4;++i){ A[i*4+i] = 1.0; }
  A[0*4+3] = dx;
  A[1*4+3] = dy;
  A[2*4+3] = dz;
}

void SetAffine_Rotate_Rodriguez
(double A[16],
 double dx, double dy, double dz)
{
  for(int i=0;i<16;++i){ A[i] = 0.0; }
  ////
  const double sqlen = dx*dx+dy*dy+dz*dz;
  const double tmp1 = 1.0/(1+0.25*sqlen);
  A[0*4+0] = 1+tmp1*(+0.5*dx*dx-0.5*sqlen);
  A[0*4+1] =  +tmp1*(-dz+0.5*dx*dy);
  A[0*4+2] =  +tmp1*(+dy+0.5*dx*dz);
  A[0*4+3] = 0.0;
  ////
  A[1*4+0] =  +tmp1*(+dz+0.5*dy*dx);
  A[1*4+1] = 1+tmp1*(+0.5*dy*dy-0.5*sqlen);
  A[1*4+2] =  +tmp1*(-dx+0.5*dy*dz);
  A[1*4+3] = 0.0;
  ////
  A[2*4+0] =  +tmp1*(-dy+0.5*dz*dx);
  A[2*4+1] =  +tmp1*(+dx+0.5*dz*dy);
  A[2*4+2] = 1+tmp1*(+0.5*dz*dz-0.5*sqlen);
  A[2*4+3] = 0.0;
  ////
  A[3*4+0] = 0.0;
  A[3*4+1] = 0.0;
  A[3*4+2] = 0.0;
  A[3*4+3] = 1.0;
}


