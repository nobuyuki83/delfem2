/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstring>
#include <cstdlib>

#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/gl.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "delfem2/gl_color.h"

static void UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3])
{
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const double invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

/*
// probably std::stroi is safer to use but it is only for C++11
static int myStoi(const std::string& str){
  char* e;
  long d = std::strtol(str.c_str(),&e,0);
  return (int)d;
}

static double myStof(const std::string& str){
  char* e;
  float fval = std::strtof(str.c_str(),&e);
  return fval;
}
 */

inline void myGlVertex3d(int i, const std::vector<double>& aV)
{
  glVertex3d(aV[i*3+0],aV[i*3+1],aV[i*3+2]);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void myGlMaterialDiffuse(const CColor& color){
  float c[4];
  c[0] = color.r;
  c[1] = color.g;
  c[2] = color.b;
  c[3] = color.a;
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

void myGlColor(const CColor& c){
  ::glColor4d(c.r, c.g, c.b, c.a );
}

void myGlColorDiffuse(const CColor& color){
  ::glColor4d(color.r, color.g, color.b, color.a );
  float c[4] = {color.r, color.g, color.b, color.a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

void myGlDiffuse(const CColor& color){
  float c[4] = {color.r, color.g, color.b, color.a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

void CColor::glColor() const {
  ::glColor4d(r, g, b, a);
}

void CColor::glMaterialDiffuse() const {
  float cf[4] = {r,g,b,a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, cf);
}

void GetRGB_HSV
(float&r, float& g, float& b,
 float h, float s, float v)
{
  r = v;
  g = v;
  b = v;
  if (s > 0.0f) {
    h *= 6.0f;
    const int i = (int) h;
    const float f = h - (float) i;
    switch (i) {
      default:
      case 0:
        g *= 1 - s * (1 - f);
        b *= 1 - s;
        break;
      case 1:
        r *= 1 - s * f;
        b *= 1 - s;
        break;
      case 2:
        r *= 1 - s;
        b *= 1 - s * (1 - f);
        break;
      case 3:
        r *= 1 - s;
        g *= 1 - s * f;
        break;
      case 4:
        r *= 1 - s * (1 - f);
        g *= 1 - s;
        break;
      case 5:
        g *= 1 - s;
        b *= 1 - s * f;
        break;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void DrawBackground(const CColor& c)
{
  glPushAttrib(GL_TRANSFORM_BIT|GL_CURRENT_BIT|GL_ENABLE_BIT);
  ::glShadeModel(GL_SMOOTH);
  GLboolean is_lighting = (glIsEnabled(GL_LIGHTING));
  GLboolean is_texture  = (glIsEnabled(GL_TEXTURE_2D));
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_TEXTURE_2D);
  ::glDisable(GL_DEPTH_TEST);
  
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glLoadIdentity();
  ::glMatrixMode(GL_PROJECTION);
  ::glPushMatrix();
  ::glLoadIdentity();
  
  ::glBegin(GL_QUADS);
  ::glColor3f(c.r,c.g,c.b);
  ::glVertex3d(-1,-1,0);
  ::glVertex3d(+1,-1,0);
  ::glColor3d(1,1,1);
  ::glVertex3d(+1,+1,0);
  ::glVertex3d(-1,+1,0);
  ::glEnd();
  
  ::glMatrixMode(GL_PROJECTION);
  ::glPopMatrix();
  ::glMatrixMode(GL_MODELVIEW);
  ::glPopMatrix();
  
  ::glPopAttrib();
  
  ::glEnable(GL_DEPTH_TEST);
  if( is_lighting ){ ::glEnable(GL_LIGHTING);   }
  if( is_texture  ){ ::glEnable(GL_TEXTURE_2D); }
}

void DrawBackground()
{
  //  ::glColor3d(0.2,0.7,0.7);
  DrawBackground( CColor(0.5, 0.5, 0.5) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void interpolateColor
(CColor& Cout, float r, const CColor& C0, const CColor& C1)
{
  Cout.r = (1-r)*C0.r+r*C1.r;
  Cout.g = (1-r)*C0.g+r*C1.g;
  Cout.b = (1-r)*C0.b+r*C1.b;
  Cout.a = (1-r)*C0.a+r*C1.a;
}

void heatmap(double input,double* color)
{
  if(0 <=input&&input <=0.25){
    color[0] = 0.0;
    color[1] = input*4.0;
    color[2] = 1.0;
  }
  else if(0.25<input && input <=0.5){
    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 2.0-input*4.0;
  }
  else if(0.5<input && input <=0.75){
    color[0] = input*4 -2.0;
    color[1] = 1.0;
    color[2] = 0.0;
  }
  else if(0.75<input&& input <=1.0){
    color[0] = 1.0;
    color[1] = 4.0 -input*4.0;
    color[2] = 0.0;
  }
  else if(1.0<input){
    color[0] = 1.0;
    color[1] = 0.0;
    color[2] = 0.0;
  }
  else{
    color[0] = 0.0;
    color[1] = 0.0;
    color[2] = 1.0;
  }
}

void heatmap_glColor(double input)
{
  double c[3]; heatmap(input,c);
  ::glColor3dv(c);
}

void heatmap_glDiffuse(double input)
{
  double c[3]; heatmap(input,c);
  float cf[4] = {(float)c[0],(float)c[1],(float)c[2],1.f};
  glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,cf);
}

CColor getColor(double input, const std::vector<std::pair<double, CColor> >& colorMap)
{
  if (colorMap.size()==0) return CColor::Black();
  if (input < colorMap[0].first){
    return colorMap[0].second;
  }
  for (int ic = 0; ic<(int)colorMap.size()-1; ++ic){
    double val0 = colorMap[ic].first;
    double val1 = colorMap[ic+1].first;
    if (val0<=input&&input<=val1){
      float rp = (float)((input-val0)/(val1-val0));
      CColor color;
      interpolateColor(color, rp, colorMap[ic].second, colorMap[ic+1].second);
      return color;
    }
  }
  return colorMap[colorMap.size()-1].second;
}

void heatmap(double input, const std::vector<std::pair<double, CColor> >& colorMap)
{
  const CColor& c = getColor(input, colorMap);
  myGlColorDiffuse(c);
}



void makeHeatMap_BlueGrayRed(std::vector<std::pair<double, CColor> >& colorMap, float min, float max)
{
  double diff = (max-min)*0.25;
  colorMap.push_back(std::make_pair(min+diff*0, CColor(0.0f, 0.0f, 1.0f, 1.0f))); // blue
  colorMap.push_back(std::make_pair(min+diff*1, CColor(0.0f, 0.2f, 1.0f, 1.0f)));
  colorMap.push_back(std::make_pair(min+diff*2, CColor(0.5f, 0.5f, 0.5f, 1.0f))); // gray
  colorMap.push_back(std::make_pair(min+diff*3, CColor(1.0f, 0.2f, 0.0f, 1.0f)));
  colorMap.push_back(std::make_pair(min+diff*4, CColor(1.0f, 0.0f, 0.0f, 1.0f))); // red
}

void makeHeatMap_BlueCyanGreenYellowRed(std::vector<std::pair<double, CColor> >& colorMap, float min, float max, float alpha)
{
  double diff = (max-min)*0.25;
  colorMap.push_back(std::make_pair(min+diff*0, CColor(0.0f, 0.0f, 1.0f, alpha))); // blue
  colorMap.push_back(std::make_pair(min+diff*1, CColor(0.0f, 1.0f, 1.0f, alpha))); // cyan
  colorMap.push_back(std::make_pair(min+diff*2, CColor(0.0f, 1.0f, 0.0f, alpha))); // green
  colorMap.push_back(std::make_pair(min+diff*3, CColor(1.0f, 1.0f, 0.0f, alpha))); // yellow
  colorMap.push_back(std::make_pair(min+diff*4, CColor(1.0f, 0.0f, 0.0f, alpha))); // red
}

void makeHeatMap_RedYellowGreenCyanBlue(std::vector<std::pair<double, CColor> >& colorMap, float min, float max)
{
  double diff = (max-min)*0.25;
  colorMap.push_back(std::make_pair(min+diff*0, CColor(1.0f, 0.0f, 0.0f, 1.0f))); // red
  colorMap.push_back(std::make_pair(min+diff*1, CColor(1.0f, 1.0f, 0.0f, 1.0f))); // yellow
  colorMap.push_back(std::make_pair(min+diff*2, CColor(0.0f, 1.0f, 0.0f, 1.0f))); // green
  colorMap.push_back(std::make_pair(min+diff*3, CColor(0.0f, 1.0f, 1.0f, 1.0f))); // cyan
  colorMap.push_back(std::make_pair(min+diff*4, CColor(0.0f, 0.0f, 1.0f, 1.0f))); // blue
}

//////////////////////////

void DrawMeshTri2D_ScalarP1
(const double* aXY, int nXY,
 const unsigned int* aTri, int nTri,
 const double* paVal,
 int nstride,
 const std::vector< std::pair<double,CColor> >& colorMap)
{
//  const unsigned int ntri = (int)aTri.size()/3;
//  const unsigned int nxys = (int)aXY.size()/2;
  glShadeModel(GL_SMOOTH);
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<nTri;++itri){
    const int ino0 = aTri[itri*3+0]; assert(ino0>=0&&ino0<nXY);
    const int ino1 = aTri[itri*3+1]; assert(ino1>=0&&ino1<nXY);
    const int ino2 = aTri[itri*3+2]; assert(ino2>=0&&ino2<nXY);
    const double v0 = paVal[ino0*nstride];
    const double v1 = paVal[ino1*nstride];
    const double v2 = paVal[ino2*nstride];
    heatmap(v0,colorMap); ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
    heatmap(v1,colorMap); ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    heatmap(v2,colorMap); ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
  }
  ::glEnd();
}

void DrawMeshTri2D_ScalarP0
(std::vector<int>& aTri,
 std::vector<double>& aXY,
 std::vector<double>& aVal,
 int nstride,
 int noffset,
 const std::vector< std::pair<double,CColor> >& colorMap)
{
  const unsigned int ntri = (int)aTri.size()/3;
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itri=0;itri<ntri;++itri){
    const int ino0 = aTri[itri*3+0];
    const int ino1 = aTri[itri*3+1];
    const int ino2 = aTri[itri*3+2];
    const double v0 = aVal[itri*nstride+noffset];
    heatmap(v0,colorMap);
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
  }
  ::glEnd();
}


void DrawSingleTri3D_Scalar_Vtx
(const double* aXYZ,
 const unsigned int* tri,
 const double* aValVtx,
 const std::vector<std::pair<double, CColor> >& colorMap)
{
  const int i0 = tri[0];
  const int i1 = tri[1];
  const int i2 = tri[2];
  if (i0==-1){
    assert(i1==-1); assert(i2==-1);
    return;
  }
//  assert(i0>=0&&i0<(int)aXYZ.size()/3);
//  assert(i1>=0&&i1<(int)aXYZ.size()/3);
//  assert(i2>=0&&i2<(int)aXYZ.size()/3);
  const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
  const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
  const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
  {
    double n[3], a; UnitNormalAreaTri3D(n, a, p0, p1, p2);
    ::glNormal3dv(n);
  }
  const double vt0 = aValVtx[i0];
  const double vt1 = aValVtx[i1];
  const double vt2 = aValVtx[i2];
  heatmap(vt0, colorMap); glVertex3dv(p0);
  heatmap(vt1, colorMap); glVertex3dv(p1);
  heatmap(vt2, colorMap); glVertex3dv(p2);
}

void DrawSingleQuad3D_Scalar_Vtx
(const std::vector<double>& aXYZ,
 const unsigned int* quad,
 const double* aValVtx,
 const std::vector<std::pair<double, CColor> >& colorMap)
{
  const int i0 = quad[0];
  const int i1 = quad[1];
  const int i2 = quad[2];
  const int i3 = quad[3];
  if (i0==-1){
    assert(i1==-1); assert(i2==-1); assert(i3==-1);
    return;
  }
  assert(i0>=0&&i0<(int)aXYZ.size()/3);
  assert(i1>=0&&i1<(int)aXYZ.size()/3);
  assert(i2>=0&&i2<(int)aXYZ.size()/3);
  assert(i3>=0&&i3<(int)aXYZ.size()/3);
  const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
  const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
  const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
  const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
  {
    double n[3], a; UnitNormalAreaTri3D(n, a, p0, p1, p2);
    ::glNormal3dv(n);
  }
  const double vt0 = aValVtx[i0];
  const double vt1 = aValVtx[i1];
  const double vt2 = aValVtx[i2];
  const double vt3 = aValVtx[i3];
  heatmap(vt0, colorMap); glVertex3dv(p0);
  heatmap(vt1, colorMap); glVertex3dv(p1);
  heatmap(vt2, colorMap); glVertex3dv(p2);
  heatmap(vt3, colorMap); glVertex3dv(p3);
}

// vetex value
void DrawMeshTri3D_ScalarP1
(const double* aXYZ, int nXYZ,
 const unsigned int* aTri, int nTri,
 const double* aValSrf,
 const std::vector<std::pair<double, CColor> >& colorMap)
{
  ::glBegin(GL_TRIANGLES);
  for (int itri = 0; itri<nTri; ++itri){
    DrawSingleTri3D_Scalar_Vtx(aXYZ, aTri+itri*3, aValSrf, colorMap);
  }
  ::glEnd();
}

// vetex value
void DrawMeshTri3D_ScalarP1
(const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 const double* aValSrf,
 const std::vector<std::pair<double, CColor> >& colorMap)
{
  const int nTri = (int)aTri.size()/3;
  const int nXYZ = (int)aXYZ.size()/3;
  DrawMeshTri3D_ScalarP1(aXYZ.data(), nXYZ,
                         aTri.data(), nTri,
                         aValSrf,
                         colorMap);
}

// vetex value
void DrawMeshElem3D_Scalar_Vtx
(const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aElemInd,
 const std::vector<unsigned int>& aElem,
 const double* aValVtx,
 const std::vector<std::pair<double, CColor> >& colorMap)
{
  if( aElemInd.size() == 0 ) return;
  /////
  for(unsigned int ielem=0;ielem<aElemInd.size()-1;++ielem){
    const int ielemind0 = aElemInd[ielem];
    const int ielemind1 = aElemInd[ielem+1];
    if( ielemind1 - ielemind0 == 3 ){
      ::glBegin(GL_TRIANGLES);
      DrawSingleTri3D_Scalar_Vtx(aXYZ.data(),
                                 aElem.data()+ielemind0,
                                 aValVtx,
                                 colorMap);
      ::glEnd();
    }
    else if(ielemind1-ielemind0 == 4){
      ::glBegin(GL_QUADS);
      DrawSingleQuad3D_Scalar_Vtx(aXYZ,
                                  aElem.data()+ielemind0,
                                  aValVtx,
                                  colorMap);
      ::glEnd();
    }
  }
}

// element-wise
void drawMeshTri3D_ScalarP0
(const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aValSrf,
 const std::vector<std::pair<double, CColor> >& colorMap)
{
  const unsigned int nTri = aTri.size()/3;
  if( aValSrf.size()!=nTri) return;
  /////
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri<nTri; ++itri){
    const int i0 = aTri[itri*3+0];
    const int i1 = aTri[itri*3+1];
    const int i2 = aTri[itri*3+2];
    if (i0==-1){
      assert(i1==-1); assert(i2==-1);
      continue;
    }
    assert(i0>=0&&i0 < (int)aXYZ.size()/3 );
    assert(i1>=0&&i1 < (int)aXYZ.size()/3 );
    assert(i2>=0&&i2 < (int)aXYZ.size()/3 );
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    double n[3], a; UnitNormalAreaTri3D(n, a, p0, p1, p2);
    const double vt = aValSrf[itri];
    ::glNormal3dv(n);
    heatmap(vt, colorMap);
    glVertex3dv(p0);
    glVertex3dv(p1);
    glVertex3dv(p2);
  }
  ::glEnd();
}



void DrawMeshTri3D_VtxColor
(const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 std::vector<CColor>& aColor)
{
  const int nTri = (int)aTri.size()/3;
  /////
  for(int itri=0;itri<nTri;++itri){
    const int i1 = aTri[itri*3+0];
    const int i2 = aTri[itri*3+1];
    const int i3 = aTri[itri*3+2];
    if( i1 == -1 ){
      assert(i2==-1); assert(i3==-1);
      continue;
    }
    ::glBegin(GL_TRIANGLES);
    assert( i1 >= 0 && i1 < (int)aXYZ.size()/3 );
    assert( i2 >= 0 && i2 < (int)aXYZ.size()/3 );
    assert( i3 >= 0 && i3 < (int)aXYZ.size()/3 );
    double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
    double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
    double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
    double un[3], area;
    UnitNormalAreaTri3D(un,area, p1,p2,p3);
    ::glNormal3dv(un);
    ::myGlColorDiffuse(aColor[i1]); ::myGlVertex3d(i1,aXYZ);
    ::myGlColorDiffuse(aColor[i2]); ::myGlVertex3d(i2,aXYZ);
    ::myGlColorDiffuse(aColor[i3]); ::myGlVertex3d(i3,aXYZ);
    ::glEnd();
  }
}

// 0: no, 1:lighting, 2:no-lighting
void DrawMeshTri3DFlag_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 const std::vector<int>& aIndGroup,
 std::vector< std::pair<int,CColor> >& aColor)
{
  const unsigned int nTri = aTri.size()/3;
  for(unsigned int itri=0;itri<nTri;++itri){
    const int ig0 = aIndGroup[itri];
    if( ig0 < 0 || ig0 >= (int)aColor.size() ) continue;
    const int imode = aColor[ig0].first;
    if(      imode == 0 ) continue;
    else if( imode == 1 ){ ::glEnable(GL_LIGHTING); }
    else if( imode == 2 ){ ::glDisable(GL_LIGHTING); }
    ::myGlColorDiffuse(aColor[ig0].second);
    const int i1 = aTri[itri*3+0];
    const int i2 = aTri[itri*3+1];
    const int i3 = aTri[itri*3+2];
    if( i1 == -1 ){
      assert(i2==-1); assert(i3==-1);
      continue;
    }
    ::glBegin(GL_TRIANGLES);
    assert( i1 >= 0 && i1 < (int)aXYZ.size()/3 );
    assert( i2 >= 0 && i2 < (int)aXYZ.size()/3 );
    assert( i3 >= 0 && i3 < (int)aXYZ.size()/3 );
    double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
    double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
    double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
    double un[3], area;
    UnitNormalAreaTri3D(un,area, p1,p2,p3);
    ::glNormal3dv(un);
    ::myGlVertex3d(i1,aXYZ);
    ::myGlVertex3d(i2,aXYZ);
    ::myGlVertex3d(i3,aXYZ);
    ::glEnd();
  }
}


//////////////////////////////////////////////////////////
// tet from here

// 3D value
void DrawMeshTet3D_ScalarP1
(const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aValSrf,
 const std::vector<std::pair<double, CColor> >& colorMap)
{
  /////
  ::glBegin(GL_TRIANGLES);
  for (int itri = 0; itri<nTet; ++itri){
    const int i0 = aTet[itri*4+0];
    const int i1 = aTet[itri*4+1];
    const int i2 = aTet[itri*4+2];
    const int i3 = aTet[itri*4+3];
    if (i0==-1){
      assert(i1==-1); assert(i2==-1);
      continue;
    }
    assert(i0>=0&&i0<nXYZ);
    assert(i1>=0&&i1<nXYZ);
    assert(i2>=0&&i2<nXYZ);
    assert(i3>=0&&i3<nXYZ);
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double un0[3], a0; UnitNormalAreaTri3D(un0,a0, p1,p2,p3);
    double un1[3], a1; UnitNormalAreaTri3D(un1,a1, p2,p0,p3);
    double un2[3], a2; UnitNormalAreaTri3D(un2,a2, p3,p0,p1);
    double un3[3], a3; UnitNormalAreaTri3D(un3,a3, p0,p2,p1);
    const double vt0 = aValSrf[i0];
    const double vt1 = aValSrf[i1];
    const double vt2 = aValSrf[i2];
    const double vt3 = aValSrf[i3];
    ::glNormal3dv(un0);
    heatmap(vt1, colorMap); glVertex3dv(p1);
    heatmap(vt2, colorMap); glVertex3dv(p2);
    heatmap(vt3, colorMap); glVertex3dv(p3);
    ::glNormal3dv(un1);
    heatmap(vt2, colorMap); glVertex3dv(p2);
    heatmap(vt3, colorMap); glVertex3dv(p3);
    heatmap(vt0, colorMap); glVertex3dv(p0);
    ::glNormal3dv(un2);
    heatmap(vt3, colorMap); glVertex3dv(p3);
    heatmap(vt0, colorMap); glVertex3dv(p0);
    heatmap(vt1, colorMap); glVertex3dv(p1);
    ::glNormal3dv(un3);
    heatmap(vt0, colorMap); glVertex3dv(p0);
    heatmap(vt1, colorMap); glVertex3dv(p1);
    heatmap(vt2, colorMap); glVertex3dv(p2);
  }
  ::glEnd();
}

static bool IsAbovePlane(const double p[3], const double org[3], const double n[3])
{
  const double dot = (p[0]-org[0])*n[0] + (p[1]-org[1])*n[1] + (p[2]-org[2])*n[2];
  return dot > 0;
}


void DrawMeshTet3D_Cut
(const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTet,
 const std::vector<CColor>& aColor,
 const double org[3], const double n[3])
{
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_DIFFUSE);
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itet=0;itet<aTet.size()/4;itet++){
    const int ino0 = aTet[itet*4+0];
    const int ino1 = aTet[itet*4+1];
    const int ino2 = aTet[itet*4+2];
    const int ino3 = aTet[itet*4+3];
    const double p0[3] = {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]};
    const double p3[3] = {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]};
    if( IsAbovePlane(p0, org, n) ) continue;
    if( IsAbovePlane(p1, org, n) ) continue;
    if( IsAbovePlane(p2, org, n) ) continue;
    if( IsAbovePlane(p3, org, n) ) continue;
    //    ::glColor3d(1,1,0);
    aColor[itet].glColorDiffuse();
    ////
    double n[3], area;
    UnitNormalAreaTri3D(n, area, p0, p2, p1);
    ::glNormal3dv(n);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    /////
    UnitNormalAreaTri3D(n, area, p0, p1, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    /////
    UnitNormalAreaTri3D(n, area, p1, p2, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
    /////
    UnitNormalAreaTri3D(n, area, p2, p0, p3);
    ::glNormal3dv(n);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  bool is_lighting = glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itet=0;itet<aTet.size()/4;itet++){
    const int ino0 = aTet[itet*4+0];
    const int ino1 = aTet[itet*4+1];
    const int ino2 = aTet[itet*4+2];
    const int ino3 = aTet[itet*4+3];
    const double p0[3] = {aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2]};
    const double p3[3] = {aXYZ[ino3*3+0], aXYZ[ino3*3+1], aXYZ[ino3*3+2]};
    if( IsAbovePlane(p0, org, n) ) continue;
    if( IsAbovePlane(p1, org, n) ) continue;
    if( IsAbovePlane(p2, org, n) ) continue;
    if( IsAbovePlane(p3, org, n) ) continue;
    ////
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p0);
    ::glVertex3dv(p2);
    ::glVertex3dv(p0);
    ::glVertex3dv(p3);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    ::glVertex3dv(p2);
    ::glVertex3dv(p3);
  }
  ::glEnd();
  ////
  /*
   ::glColor3d(0,0,0);
   ::glPointSize(3);
   ::glBegin(GL_POINTS);
   for(unsigned int ino=0;ino<nXYZ_;ino++){
   ::glVertex3dv(pXYZ_+ino*3);
   }
   ::glEnd();
   */
  if( is_lighting ){ glEnable(GL_LIGHTING); }
}



/////////////

void Write_Ply_Tri2DMesh_HeightColor
(const std::string& fname,
 const std::vector<int>& aTri1,
 const std::vector<double>& aXY1,
 const std::vector<double>& aVal,
 std::vector< std::pair<double,CColor> >& colorMap)
{
  const int np = aXY1.size()/2;
  const int ntri = aTri1.size()/3;
  std::ofstream fout;
  fout.open(fname.c_str(),std::ios::out);
  fout << "ply" << std::endl;
  fout << "format ascii 1.0" << std::endl;
  fout << "element vertex " << np << std::endl;
  fout << "property float x" << std::endl;
  fout << "property float y" << std::endl;
  fout << "property float z" << std::endl;
  fout << "property uchar red" << std::endl;
  fout << "property uchar green" << std::endl;
  fout << "property uchar blue" << std::endl;
  fout << "element face " << ntri << std::endl;
  fout << "property list uchar int vertex_indices" << std::endl;
  fout << "end_header" << std::endl;
  for(int ip=0;ip<np;++ip){
    double v = aVal[ip];
    CColor c = getColor(v,colorMap);
    int cr, cg, cb;
    c.getRGBChar(cr,cg,cb);
    fout << aXY1[ip*2+0] << " " << aXY1[ip*2+1] << " " << v << " ";
    fout << cr << " " << cg << " " << cb << std::endl;
  }
  for(int itri=0;itri<ntri;++itri){
    fout << "3 " << aTri1[itri*3+0] << " " << aTri1[itri*3+1] << " " << aTri1[itri*3+2] << std::endl;
  }
}

/////////////

