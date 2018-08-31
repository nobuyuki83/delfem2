#ifndef COLOR_GL_H
#define COLOR_GL_H

// OpenGL dependency (not GLUT)
// this file should not depend anything other than OpenGL

#if defined(__APPLE__) && defined(__MACH__)
  #include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
  #include <GL/gl.h>
#elif defined(WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else
  #include <GL/gl.h>
#endif

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
  void setRandomColor(){
    r = (float)rand()/(RAND_MAX+1.0f);
    g = (float)rand()/(RAND_MAX+1.0f);
    b = (float)rand()/(RAND_MAX+1.0f);
  }
  void setRandomVividColor(){
    double h = (double)rand()/(RAND_MAX+1.0);
    GetRGB_HSV(r,g,b, (float)h,1,1);
  }
  void glColor() const {
    ::glColor4d(r, g, b, a);
  }
  void glMaterialDiffuse() const {
    float cf[4] = {r,g,b,a};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, cf);
  }
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

inline void myGlMaterialDiffuse(const CColor& color){
  float c[4];
  c[0] = color.r;
  c[1] = color.g;
  c[2] = color.b;
  c[3] = color.a;
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

inline void myGlColor(const CColor& c){
  ::glColor4d(c.r, c.g, c.b, c.a );
}

inline void myGlColorDiffuse(const CColor& color){
  ::glColor4d(color.r, color.g, color.b, color.a );
  float c[4] = {color.r, color.g, color.b, color.a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

inline void myGlDiffuse(const CColor& color){
  float c[4] = {color.r, color.g, color.b, color.a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

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
                          const std::vector<int>& aTri,
                          const std::vector<double>& aValSrf,
                          const std::vector<std::pair<double, CColor> >& colorMap);

void DrawMeshTri2D_ScalarP1(std::vector<int>& aTri,
                        std::vector<double>& aXY,
                        std::vector<double>& aVal,
                        int nstride,
                        int noffset,
                        const std::vector< std::pair<double,CColor> >& colorMap);

void DrawMeshTri2D_ScalarP0(std::vector<int>& aTri,
                        std::vector<double>& aXY,
                        std::vector<double>& aVal,
                        int nstride,
                        int noffset,
                        const std::vector< std::pair<double,CColor> >& colorMap);

// 3D value -- vtx value
void DrawMeshTri3D_ScalarP1(const std::vector<double>& aXYZ,
                        const std::vector<int>& aTri,
                        const std::vector<double>& aValSrf,
                        const std::vector<std::pair<double, CColor> >& colorMap);

// scalar value on 3D mesh (mixed elem).
void DrawMeshElem3D_Scalar_Vtx(const std::vector<double>& aXYZ,
                               const std::vector<int>& aElemInd,
                               const std::vector<int>& aElem,
                               const std::vector<double>& aValVtx,
                               const std::vector<std::pair<double, CColor> >& colorMap);

// 3D value
void DrawMeshTet3D_ScalarP1(const std::vector<double>& aXYZ,
                        const std::vector<int>& aTet,
                        const std::vector<double>& aValSrf,
                        const std::vector<std::pair<double, CColor> >& colorMap);

////////////////////////////////////

void Write_Ply_Tri2DMesh_HeightColor(const std::string& fname,
                                     const std::vector<int>& aTri1,
                                     const std::vector<double>& aXY1,
                                     const std::vector<double>& aVal,
                                     std::vector< std::pair<double,CColor> >& colorMap);


////////////////////////////////////


class CMaterial{
public:
  std::string name_mtl;
  float Kd[4];
  float Ka[4];
  float Ks[4];
  float Ke[4];
  float Ns;
  int illum;
  std::string map_Kd;
public:
  void GL() const;
};

void Load_Mtl(const std::string& fname,
              std::vector<CMaterial>& aMtl);

#endif
