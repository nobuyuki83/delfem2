/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef COLOR_GL_H
#define COLOR_GL_H

// OpenGL dependency (not GLUT)
// this file should not depend anything other than OpenGL

#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <assert.h>

void GetRGB_HSV(float&r, float& g, float& b,
                float h, float s, float v);

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
class CColor
{
public:
  CColor(){
    this->r = 0.8f;
    this->g = 0.8f;
    this->b = 0.8f;
    this->a = 1.0f;
  }
  CColor(float r, float g, float b){
    this->r = r;
    this->g = g;
    this->b = b;
    a = 1.0;
  }
  CColor(float r, float g, float b, float a){
    this->r = r;
    this->g = g;
    this->b = b;
    this->a = a;
  }
  CColor(const std::vector<double>& v){
    this->r = 1.0;
    this->g = 1.0;
    this->b = 1.0;
    this->a = 1.0;
    if( v.size() == 4 ){
      this->r = v[0];
      this->g = v[1];
      this->b = v[2];
      this->a = v[3];
    }
    else if( v.size() == 3 ){
      this->r = v[0];
      this->g = v[1];
      this->b = v[2];
    }
    else if( v.size() > 0 ){
      this->r = v[0];
      this->g = v[0];
      this->b = v[0];
    }
  }
  void setRandomColor(){
    r = (float)rand()/(RAND_MAX+1.0f);
    g = (float)rand()/(RAND_MAX+1.0f);
    b = (float)rand()/(RAND_MAX+1.0f);
  }
  void setRandomVividColor(){
    double h = (double)rand()/(RAND_MAX+1.0);
    GetRGB_HSV(r,g,b, (float)h,1,1);
  }
  void glColor() const;
  void glMaterialDiffuse() const;
  void glColorDiffuse() const {
    this->glColor();
    this->glMaterialDiffuse();
  }
  void getRGBChar(int& cr, int& cg, int& cb) const {
    if( r < 0 ){ cr = 0; }
    else if( r > 1 ){ cr = 255; }
    else{ cr = (int)(255*r); }
    ///
    if( g < 0 ){ cg = 0; }
    else if( g > 1 ){ cg = 255; }
    else{ cg = (int)(255*g); }
    ////
    if( b < 0 ){ cb = 0; }
    else if( b > 1 ){ cb = 255; }
    else{ cb = (int)(255*b); }
  }
  static CColor Black() { return CColor(0,0,0); }
  static CColor Red() { return CColor(1,0,0); }
  static CColor Blue() { return CColor(0,0,1); }
  static CColor Green(){ return CColor(0,1,0); }
  static CColor Yellow(){ return CColor(1,1,0); }
  static CColor Orange(){ return CColor(1, 0.71f, 0.30f); }
  static CColor Purple(){ return CColor(1,0,1); }
  static CColor Cyan(){ return CColor(0,1,1); }
  static CColor White() { return CColor(1,1,1); }
  static CColor Gray(float f) { return CColor(f, f, f); }
  static CColor Gray() { return CColor(0.8f, 0.8f, 0.8f); }
public:
  float r;
  float g;
  float b;
  float a;
};

void interpolateColor(CColor& Cout, float r, const CColor& C0, const CColor& C1);

void myGlMaterialDiffuse(const CColor& color);
void myGlColor(const CColor& c);
void myGlColorDiffuse(const CColor& color);
inline void myGlDiffuse(const CColor& color);

void DrawBackground(const CColor& c);
void DrawBackground();

void heatmap(double input,double* color);
void heatmap_glColor(double input);
void heatmap_glDiffuse(double input);

////////////////////////////////////////////////////////////////////////////////////
// vector from here

void heatmap(double input, const std::vector<std::pair<double, CColor> >& colorMap);
CColor getColor(double input, const std::vector<std::pair<double, CColor> >& colorMap);


void makeHeatMap_BlueGrayRed(std::vector<std::pair<double, CColor> >& colorMap, float min, float max);
void makeHeatMap_BlueCyanGreenYellowRed(std::vector<std::pair<double, CColor> >& colorMap, float min, float max, float alpha=1);
void makeHeatMap_RedYellowGreenCyanBlue(std::vector<std::pair<double, CColor> >& colorMap, float min, float max);

class CColorMap
{
public:
  CColorMap(){}
  CColorMap(double min, double max, const std::string& str){
    if( str == "bgr" ){
      makeHeatMap_BlueGrayRed(aColor, min, max);
    }
    else{
      makeHeatMap_BlueCyanGreenYellowRed(aColor, min, max);
    }
  }
public:
  std::vector< std::pair<double,CColor> > aColor;
};

// 0: no, 1:lighting, 2:no-lighting
void DrawMeshTri3D_VtxColor(const std::vector<double>& aXYZ,
                            const std::vector<int>& aTri,
                            std::vector<CColor>& aColor);

// 0: no, 1:lighting, 2:no-lighting
void DrawMeshTri3DFlag_FaceNorm(const std::vector<double>& aXYZ,
                            const std::vector<int>& aTri,
                            const std::vector<int>& aIndGroup,
                            std::vector< std::pair<int,CColor> >& aColor);

void DrawMeshTri_ScalarP0(const std::vector<double>& aXYZ,
                          const std::vector<unsigned int>& aTri,
                          const std::vector<double>& aValSrf,
                          const std::vector<std::pair<double, CColor> >& colorMap);

void DrawMeshTri2D_ScalarP1(const double* aXY, int nXY,
                            const unsigned int* aTri, int nTri,
                            const double* aVal,
                            int nstride,
                          const std::vector< std::pair<double,CColor> >& colorMap);

void DrawMeshTri2D_ScalarP0(std::vector<int>& aTri,
                            std::vector<double>& aXY,
                            std::vector<double>& aVal,
                            int nstride,
                            int noffset,
                            const std::vector< std::pair<double,CColor> >& colorMap);

// 3D value -- vtx value
void DrawMeshTri3D_ScalarP1(const double* aXYZ, int nXYZ,
                            const unsigned int* aTri, int nTri,
                            const double* aValSrf,
                            const std::vector<std::pair<double, CColor> >& colorMap);
void DrawMeshTri3D_ScalarP1(const std::vector<double>& aXYZ,
                        const std::vector<unsigned int>& aTri,
                        const double* aValSrf,
                        const std::vector<std::pair<double, CColor> >& colorMap);

// scalar value on 3D mesh (mixed elem).
void DrawMeshElem3D_Scalar_Vtx(const std::vector<double>& aXYZ,
                               const std::vector<unsigned int>& aElemInd,
                               const std::vector<unsigned int>& aElem,
                               const double* aValVtx,
                               const std::vector<std::pair<double, CColor> >& colorMap);

// 3D value
void DrawMeshTet3D_ScalarP1(const double* aXYZ, int nXYZ,
                            const unsigned int* aTet, int nTet,
                            const double* aValSrf,
                            const std::vector<std::pair<double, CColor> >& colorMap);

void DrawMeshTet3D_Cut(const std::vector<double>& aXYZ,
                       const std::vector<unsigned int>& aTet,
                       const std::vector<CColor>& aColor,
                       const double org[3], const double n[3]);

////////////////////////////////////

void Write_Ply_Tri2DMesh_HeightColor(const std::string& fname,
                                     const std::vector<int>& aTri1,
                                     const std::vector<double>& aXY1,
                                     const std::vector<double>& aVal,
                                     std::vector< std::pair<double,CColor> >& colorMap);



#endif
