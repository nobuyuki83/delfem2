/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef CAMERA_GL_H
#define CAMERA_GL_H

#include <math.h>
#include <stdio.h> // memcpy
#include <string.h>

#if defined(__APPLE__) && defined(__MACH__)
  #include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
  #include <GL/glu.h>
#elif defined(WIN32)
  #include <windows.h>
  #include <GL/glu.h>
#else
  #include <GL/glu.h>
#endif

/////////////////////////////////////////////////////////////

void getPosOnScreen_Camera2D(double& x, double& y,
                             int i, int j);
void setGL_Camera2D();

/////////////////////////////////////////////////////////////

void glhFrustumf2(float *matrix, float left, float right, float bottom, float top,
                  float znear, float zfar);
void glhPerspectivef2(float *matrix, float fovyInDegrees, float aspectRatio,
                      float znear, float zfar);
void glhTranslatef2(float *matrix, float x, float y, float z);
void glhLookAtf2(float *matrix,
                 float eyex, float eyey, float eyez,
                 float cntx, float cnty, float cntz,
                 float upx, float upy, float upz );

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

enum CAMERA_ROT_MODE
{
  CAMERA_ROT_YTOP,
  CAMERA_ROT_ZTOP,
  CAMERA_ROT_TBALL
};

class CCamera
{
public:
  CCamera(){
    is_pars = false;
    fovy = 60;
    
    view_height = 0.1;
    scale = 1.0;
    
    trans[0] = 0;
    trans[1] = 0;
    trans[2] = 0;
    
    camera_rot_mode = CAMERA_ROT_YTOP;

    psi = 0;
    theta = 0;
    
    Quat_tball[0]=1;
    Quat_tball[1]=0;
    Quat_tball[2]=0;
    Quat_tball[3]=0;
  }
  void SetGL_Camera(int win_w, int win_h);
  void Scale(double s);
  void Rot_Camera(double dx, double dy);
  void Pan_Camera(double dx, double dy);
private:
public:
  bool is_pars;
  double fovy;
  double scale;
  double view_height;
  double trans[3];
  
  CAMERA_ROT_MODE camera_rot_mode;
  // ytop
  double theta;
  double psi;
  // tball
  double Quat_tball[4];
};


std::ostream &operator<<(std::ostream &output, CCamera& c);
std::istream &operator>>(std::istream &input, CCamera& c);



#endif
