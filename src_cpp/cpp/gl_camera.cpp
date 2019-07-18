/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cmath>
#include <cstdio> // memcpy
#include <cstring>
#include <iostream>

#if defined(__APPLE__) && defined(__MACH__)
  #include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
  #include <GL/glu.h>
#elif defined(_WIN32)
  #include <windows.h>
  #include <GL/glu.h>
#else
  #include <GL/glu.h>
#endif

#include "delfem2/gl_camera.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// multiply two quaternion
static void QuatMult(double r[], const double p[], const double q[])
{
  r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
  r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
  r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
  r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}
/*
static void QuatSetIdentity(double q[]){
  q[0] = 1;
  q[1] = 0;
  q[2] = 0;
  q[3] = 0;
}
 */

// transform vector with quaternion
inline void QuatVec(double vo[], const double q[], const double vi[])
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
  
  vo[0] = (1.0 - y2 - z2)*vi[0] + (xy + zw      )*vi[1] + (zx - yw      )*vi[2];
  vo[1] = (xy - zw      )*vi[0] + (1.0 - z2 - x2)*vi[1] + (yz + xw      )*vi[2];
  vo[2] = (zx + yw      )*vi[0] + (yz - xw      )*vi[1] + (1.0 - x2 - y2)*vi[2];
}

// transform vector with conjugate of quaternion
inline void QuatConjVec(double vo[], const double q[], const double vi[])
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
  
  vo[0] = (1.0 - y2 - z2)*vi[0] + (xy - zw      )*vi[1] + (zx + yw      )*vi[2];
  vo[1] = (xy + zw      )*vi[0] + (1.0 - z2 - x2)*vi[1] + (yz - xw      )*vi[2];
  vo[2] = (zx - yw      )*vi[0] + (yz + xw      )*vi[1] + (1.0 - x2 - y2)*vi[2];
}

// copy quaternion
static void CopyQuat(double r[], const double p[])
{
  r[0] = p[0];
  r[1] = p[1];
  r[2] = p[2];
  r[3] = p[3];
}

static void RotTransAffineMatQuaternion(double r[], const double q[])
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
  
  r[ 0] = 1.0 - y2 - z2;
  r[ 1] = xy + zw;
  r[ 2] = zx - yw;
  r[ 4] = xy - zw;
  r[ 5] = 1.0 - z2 - x2;
  r[ 6] = yz + xw;
  r[ 8] = zx + yw;
  r[ 9] = yz - xw;
  r[10] = 1.0 - x2 - y2;
  r[ 3] = r[ 7] = r[11] = r[12] = r[13] = r[14] = 0.0;
  r[15] = 1.0;
}


/////////////////////////////////////////////////////////////

void getPosOnScreen_Camera2D
(double& x, double& y,
 int i, int j)
{
  int viewport[8];
  glGetIntegerv(GL_VIEWPORT, viewport);
  double hw = (double)viewport[2]*0.5; // half width
  double hh = (double)viewport[3]*0.5; // half height
  double asp = hw/hh;
  x = (i-hw)/hw*asp;
  y = (hh-j)/hh;
}

void setGL_Camera2D()
{
  int viewport[8];
  glGetIntegerv(GL_VIEWPORT, viewport);
  double w = (double)viewport[2];
  double h = (double)viewport[3];
  double asp = w/h;
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::glOrtho(-asp*2, +asp*2, -2, +2, -10, +10);
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
}

/////////////////////////////////////////////////////////////

void glhFrustumf2
(float *matrix,
 float left, float right,
 float bottom, float top,
 float znear, float zfar)
{
  float temp, temp2, temp3, temp4;
  temp = 2.0 * znear;
  temp2 = right - left;
  temp3 = top - bottom;
  temp4 = zfar - znear;
  matrix[0] = temp / temp2;
  matrix[1] = 0.0;
  matrix[2] = 0.0;
  matrix[3] = 0.0;
  matrix[4] = 0.0;
  matrix[5] = temp / temp3;
  matrix[6] = 0.0;
  matrix[7] = 0.0;
  matrix[8] = (right + left) / temp2;
  matrix[9] = (top + bottom) / temp3;
  matrix[10] = (-zfar - znear) / temp4;
  matrix[11] = -1.0;
  matrix[12] = 0.0;
  matrix[13] = 0.0;
  matrix[14] = (-temp * zfar) / temp4;
  matrix[15] = 0.0;
}

//matrix will receive the calculated perspective matrix.
//You would have to upload to your shader
// or use glLoadMatrixf if you aren't using shaders.
void glhPerspectivef2
(float *matrix,
 float fovyInDegrees, float aspectRatio,
 float znear, float zfar)
{
  float ymax, xmax;
//  float temp, temp2, temp3, temp4;
  ymax = znear * tanf(fovyInDegrees * 3.14159 / 360.0);
  //ymin = -ymax;
  //xmin = -ymax * aspectRatio;
  xmax = ymax * aspectRatio;
  glhFrustumf2(matrix, -xmax, xmax, -ymax, ymax, znear, zfar);
}


static void Normalize3D(float vec[3])
{
  float len = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  float leninv = 1.0/len;
  vec[0] *= leninv;
  vec[1] *= leninv;
  vec[2] *= leninv;
}

static void Cross3D(float r[3], const float v1[3], const float v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

void glhTranslatef2
(float *matrix,
 float x, float y, float z)
{
  matrix[12]=matrix[0]*x+matrix[4]*y+matrix[8]*z+matrix[12];
  matrix[13]=matrix[1]*x+matrix[5]*y+matrix[9]*z+matrix[13];
  matrix[14]=matrix[2]*x+matrix[6]*y+matrix[10]*z+matrix[14];
  matrix[15]=matrix[3]*x+matrix[7]*y+matrix[11]*z+matrix[15];
}

//PURPOSE:      For square matrices. This is column major for OpenGL
inline void MultiplyMatrices4by4OpenGL_FLOAT(float *result, float *matrix1, float *matrix2)
{
  result[0]=matrix1[0]*matrix2[0]+
  matrix1[4]*matrix2[1]+
  matrix1[8]*matrix2[2]+
  matrix1[12]*matrix2[3];
  result[4]=matrix1[0]*matrix2[4]+
  matrix1[4]*matrix2[5]+
  matrix1[8]*matrix2[6]+
  matrix1[12]*matrix2[7];
  result[8]=matrix1[0]*matrix2[8]+
  matrix1[4]*matrix2[9]+
  matrix1[8]*matrix2[10]+
  matrix1[12]*matrix2[11];
  result[12]=matrix1[0]*matrix2[12]+
  matrix1[4]*matrix2[13]+
  matrix1[8]*matrix2[14]+
  matrix1[12]*matrix2[15];
  
  result[1]=matrix1[1]*matrix2[0]+
  matrix1[5]*matrix2[1]+
  matrix1[9]*matrix2[2]+
  matrix1[13]*matrix2[3];
  result[5]=matrix1[1]*matrix2[4]+
  matrix1[5]*matrix2[5]+
  matrix1[9]*matrix2[6]+
  matrix1[13]*matrix2[7];
  result[9]=matrix1[1]*matrix2[8]+
  matrix1[5]*matrix2[9]+
  matrix1[9]*matrix2[10]+
  matrix1[13]*matrix2[11];
  result[13]=matrix1[1]*matrix2[12]+
  matrix1[5]*matrix2[13]+
  matrix1[9]*matrix2[14]+
  matrix1[13]*matrix2[15];
  
  result[2]=matrix1[2]*matrix2[0]+
  matrix1[6]*matrix2[1]+
  matrix1[10]*matrix2[2]+
  matrix1[14]*matrix2[3];
  result[6]=matrix1[2]*matrix2[4]+
  matrix1[6]*matrix2[5]+
  matrix1[10]*matrix2[6]+
  matrix1[14]*matrix2[7];
  result[10]=matrix1[2]*matrix2[8]+
  matrix1[6]*matrix2[9]+
  matrix1[10]*matrix2[10]+
  matrix1[14]*matrix2[11];
  result[14]=matrix1[2]*matrix2[12]+
  matrix1[6]*matrix2[13]+
  matrix1[10]*matrix2[14]+
  matrix1[14]*matrix2[15];
  
  result[3]=matrix1[3]*matrix2[0]+
  matrix1[7]*matrix2[1]+
  matrix1[11]*matrix2[2]+
  matrix1[15]*matrix2[3];
  result[7]=matrix1[3]*matrix2[4]+
  matrix1[7]*matrix2[5]+
  matrix1[11]*matrix2[6]+
  matrix1[15]*matrix2[7];
  result[11]=matrix1[3]*matrix2[8]+
  matrix1[7]*matrix2[9]+
  matrix1[11]*matrix2[10]+
  matrix1[15]*matrix2[11];
  result[15]=matrix1[3]*matrix2[12]+
  matrix1[7]*matrix2[13]+
  matrix1[11]*matrix2[14]+
  matrix1[15]*matrix2[15];
}

void glhLookAtf2
( float *matrix,
 float eyex, float eyey, float eyez,
 float cntx, float cnty, float cntz,
 float upx, float upy, float upz )
{
  float eyePosition3D[3] = {eyex, eyey, eyez};
  float center3D[3] = {cntx, cnty, cntz};
  float upVector3D[3] = {upx, upy, upz};
  
  
  //------------------
  float forward[3] = {0,0,1};
  forward[0] = center3D[0] - eyePosition3D[0];
  forward[1] = center3D[1] - eyePosition3D[1];
  forward[2] = center3D[2] - eyePosition3D[2];
  Normalize3D(forward);
  //------------------
  //Side = forward x up
  float side[3] = {1,0,0};
  Cross3D(side, forward, upVector3D);
  Normalize3D(side);
  //------------------
  //Recompute up as: up = side x forward
  float up[3] = {0,1,0};
  Cross3D(up, side, forward);
  //------------------
  float matrix2[16];
  matrix2[0] = side[0];
  matrix2[4] = side[1];
  matrix2[8] = side[2];
  matrix2[12] = 0.0;
  //------------------
  matrix2[1] = up[0];
  matrix2[5] = up[1];
  matrix2[9] = up[2];
  matrix2[13] = 0.0;
  //------------------
  matrix2[2] = -forward[0];
  matrix2[6] = -forward[1];
  matrix2[10] = -forward[2];
  matrix2[14] = 0.0;
  //------------------
  matrix2[3] = matrix2[7] = matrix2[11] = 0.0;
  matrix2[15] = 1.0;
  //------------------
  float resultMatrix[16];
  MultiplyMatrices4by4OpenGL_FLOAT(resultMatrix, matrix, matrix2);
  glhTranslatef2(resultMatrix,
                 -eyePosition3D[0], -eyePosition3D[1], -eyePosition3D[2]);
  //------------------
  memcpy(matrix, resultMatrix, 16*sizeof(float));
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void CCamera::SetGL_Camera(int win_w, int win_h)
{
    double depth = view_height/(scale*tan(0.5*fovy*3.1415/180.0));
    {
      double asp = (double)win_w/win_h;
      ::glMatrixMode(GL_PROJECTION);
      ::glLoadIdentity();
      if( is_pars ){
//        ::gluPerspective(fovy, asp, depth*0.01, depth*10);
        float mP[16];
        glhPerspectivef2(mP, fovy, asp, depth*0.01, depth*10);
        ::glMultMatrixf(mP);
      }
      else{
        ////
        ::glOrtho(-view_height/scale*asp,
                  +view_height/scale*asp,
                  -view_height/scale,
                  +view_height/scale,
                  -depth*10,
                  +depth*10);
//                  -view_height/scale*10.0,
//                  +view_height/scale*10.0);
      }
    }
    //// solve translation rotation from here
    {
      ::glMatrixMode(GL_MODELVIEW);
      ::glLoadIdentity();
      ::glTranslated(trans[0],trans[1],-depth);
    }
    if(      camera_rot_mode == CAMERA_ROT_YTOP  ){
      double x = sin(theta);
      double z = cos(theta);
      double y = sin(psi);
      x *= cos(psi);
      z *= cos(psi);
      float mMV[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
      glhLookAtf2(mMV, x,y,z, 0,0,0, 0,1,0);
      ::glMultMatrixf(mMV);
//      ::gluLookAt(x,y,z, 0,0,0, 0,1,0);
    }
    else if( camera_rot_mode == CAMERA_ROT_ZTOP  ){
      double x = sin(theta);
      double y = cos(theta);
      double z = sin(psi);
      x *= cos(psi);
      y *= cos(psi);
      float mMV[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
      glhLookAtf2(mMV, x,y,z, 0,0,0, 0,0,1);
      ::glMultMatrixf(mMV);
//      ::gluLookAt(x,y,z, 0,0,0, 0,0,1);
    }
    else if( camera_rot_mode == CAMERA_ROT_TBALL ){
      double Rview[16];
      RotTransAffineMatQuaternion(Rview,Quat_tball);
      ::glMultMatrixd(Rview);
    }
}

void CCamera::Scale(double s){
  scale *= s;
}
void CCamera::Rot_Camera(double dx, double dy){
    if(      camera_rot_mode == CAMERA_ROT_YTOP ){
      theta -= dx;
      psi   -= dy;
    }
    else if( camera_rot_mode == CAMERA_ROT_ZTOP ){
      theta += dx;
      psi   -= dy;
    }
    else if( camera_rot_mode == CAMERA_ROT_TBALL ){
      double a = sqrt(dx * dx + dy * dy);
      double ar = a*0.5; // angle
      double dq[4] = { cos(ar), -dy*sin(ar)/a, dx*sin(ar)/a, 0.0 };
      if (a != 0.0) {
        double qtmp[4]; QuatMult(qtmp, dq, Quat_tball);
        CopyQuat(Quat_tball,qtmp);
      }
    }
}

void CCamera::Pan_Camera(double dx, double dy){
    double s = view_height/scale;
    trans[0] += s*dx;
    trans[1] += s*dy;
    trans[2] += s*0.0;
  
}

std::ostream &operator<<(std::ostream &output, CCamera& c)
{
  output.setf(std::ios::scientific);
  output << c.is_pars << std::endl;
  output << c.fovy << std::endl;
  output << c.scale << std::endl;
  output << c.view_height << std::endl;
  output << c.trans[0] << " " << c.trans[1] << " " << c.trans[2] << std::endl;
  output << (int)c.camera_rot_mode << std::endl;
  output << c.theta << " " << c.psi << std::endl;
  output << c.Quat_tball[0] << " " << c.Quat_tball[1] << " " << c.Quat_tball[2] << " " << c.Quat_tball[3] << std::endl;
  return output;
}

std::istream &operator>>(std::istream &input, CCamera& c)
{
  {
    int is0; input >> is0; c.is_pars = (bool)is0;
  }
  input >> c.fovy;
  input >> c.scale;
  input >> c.view_height;
  input >> c.trans[0] >> c.trans[1] >> c.trans[2];
  {
    int imode; input >> imode; c.camera_rot_mode = (CAMERA_ROT_MODE)imode;
  }
  input >> c.theta >> c.psi;
  input >> c.Quat_tball[0] >> c.Quat_tball[1] >> c.Quat_tball[2] >> c.Quat_tball[3];
  return input;
}

