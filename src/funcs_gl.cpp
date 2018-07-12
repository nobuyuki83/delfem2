#include <stdio.h>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#if defined(__APPLE__) && defined(__MACH__) // Mac
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
  #include <GL/glu.h>
#elif defined(WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
  #include <GL/glu.h>
#else // linux
  #include <GL/gl.h>
  #include <GL/glu.h>
#endif

#include "delfem2/funcs_gl.h"

void DrawAxis(double s)
{
  GLboolean is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(1,0,0);
  ::glVertex3d(0,0,0);
  ::glVertex3d(s,0,0);
  ::glColor3d(0,1,0);
  ::glVertex3d(0,0,0);
  ::glVertex3d(0,s,0);
  ::glColor3d(0,0,1);
  ::glVertex3d(0,0,0);
  ::glVertex3d(0,0,s);
  ::glEnd();
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
}


void DrawSphere
(int nla, int nlo)
{
  if( nla <= 1 || nlo <= 2 ){ return; }
  const double pi = 3.1415926535;
  double dla = 2.0*pi/nla;
  double dlo = pi/nlo;
  ::glBegin(GL_QUADS);
  for(int ila=0;ila<nla;ila++){
    for(int ilo=0;ilo<nlo;ilo++){
      double rla0 = (ila+0)*dla;
      double rla1 = (ila+1)*dla;
      double rlo0 = (ilo+0)*dlo;
      double rlo1 = (ilo+1)*dlo;
      double p0[3] = { cos(rla0)*cos(rlo0), cos(rla0)*sin(rlo0), sin(rla0) };
      double p1[3] = { cos(rla0)*cos(rlo1), cos(rla0)*sin(rlo1), sin(rla0) };
      double p2[3] = { cos(rla1)*cos(rlo1), cos(rla1)*sin(rlo1), sin(rla1) };
      double p3[3] = { cos(rla1)*cos(rlo0), cos(rla1)*sin(rlo0), sin(rla1) };
      ::glVertex3dv(p0);
      ::glVertex3dv(p1);
      ::glVertex3dv(p2);
      ::glVertex3dv(p3);
    }
  }
  ::glEnd();
}

void DrawSphere(int nla, int nlo, double rad, double x, double y, double z)
{
  ::glTranslated(+x,+y,+z);
  ::glScaled(rad, rad, rad);
  DrawSphere(nla,nlo);
  ::glScaled(1.0/rad, 1.0/rad, 1.0/rad);
  ::glTranslated(-x,-y,-z);
}

void setSomeLighting()
{
  glEnable(GL_LIGHTING);
  //  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, 1.0);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 0.0);
//  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  {
    glEnable(GL_LIGHT0);
    GLfloat light0_Kd[]  = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat light0_Specular[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat light0_Pos[4] = {0.25f, 1.0f, +1.25f, 0.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_Kd);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_Specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_Pos);
  }
  {
    glEnable(GL_LIGHT1);
    GLfloat light1_Kd[]  = {0.3f, 0.3f, 0.3f, 1.0f};
    GLfloat light1_Pos[4] = {0.00f, 0.0f, +1.00f, 0.0f};
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_Kd);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_Pos);
  }
}

void setSomeLighting2()
{
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
  { // initialize light parameter
    GLfloat light0_Kd[]   = {0.9f, 0.3f, 0.3f, 1.0f};
    GLfloat light0_Pos[4] = {+0.5f, -0.5f, +1.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0_Kd);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_Pos);
    ////
    GLfloat light1_Kd[]   = {0.3f, 0.3f, 0.9f, 1.0f};
    GLfloat light1_Pos[4] = {-0.5f, +0.5f, +1.0f, 0.0f};
    glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_Kd);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_Pos);
  }
}

void setSomeLighting3()
{
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 0.0);
  {
    glEnable(GL_LIGHT0);
    GLfloat light0_Kd[]  = {0.8f, 0.8f, 0.8f, 1.0f};
    GLfloat light0_Specular[4] = {0.5f, 0.5f, .5f, 1.0f};
    GLfloat light0_Pos[4] = {0.25f, 1.0f, +1.25f, 0.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_Kd);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_Specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_Pos);
  }
  {
    glEnable(GL_LIGHT1);
    GLfloat light1_Kd[]  = {0.3f, 0.3f, 0.3f, 1.0f};
    GLfloat light1_Pos[4] = {-1.00f, 0.0f, +1.00f, 0.0f};
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_Kd);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_Pos);
  }
}

static void ShadowMatrix(float m[16], float plane[4], float lpos[3])
{
  float dot = plane[0]*lpos[0] + plane[1]*lpos[1] + plane[2]*lpos[2] + plane[3];
  for(int j=0; j<4;++j){
    for(int i=0; i<4; ++i){
      m[j*4+i] = - plane[j]*lpos[i];
      if (i == j){ m[j*4+i] += dot; }
    }
  }
}

void drawFloorShadow(void (*DrawObject)(), float yfloor, float wfloor)
{
  GLboolean is_lighting = ::glIsEnabled(GL_LIGHTING);
  {
    ::glClearStencil(0);
    { // draw floor (stencil 1)
      glEnable(GL_STENCIL_TEST);
      glStencilFunc( GL_ALWAYS, 1, ~0);
      glStencilOp(GL_KEEP,GL_KEEP ,GL_REPLACE);
      { // floor
        ::glDisable(GL_LIGHTING);
        glColor4f(0.6f, 0.6f, 0.5f, 1.0f);
        ::glBegin(GL_QUADS);
        ::glNormal3d(0,1,0);
        ::glVertex3d(-wfloor,yfloor,-wfloor);
        ::glVertex3d(+wfloor,yfloor,-wfloor);
        ::glVertex3d(+wfloor,yfloor,+wfloor);
        ::glVertex3d(-wfloor,yfloor,+wfloor);
        ::glEnd();
      }
    }
    { // draw stensil
      glColorMask(0,0,0,0);
      glDepthMask(0);
      glEnable(GL_STENCIL_TEST);
      glStencilFunc( GL_EQUAL, 1, ~0);
      glStencilOp(GL_KEEP,GL_KEEP ,GL_INCR);
      glPushMatrix();
      {
        float plane[4] = {0,1,0,-yfloor-0.001f};
        float lpos[4] = {0,5,0,1};
        float m_shadow[16]; ShadowMatrix(m_shadow, plane, lpos);
        glMultMatrixf(m_shadow);
      }
      DrawObject();
      glPopMatrix();
      glColorMask(1,1,1,1);
      glDepthMask(1);
    }
    { // draw shadow
      glStencilFunc( GL_EQUAL, 2, ~0 );
      glStencilOp(GL_KEEP, GL_KEEP ,GL_KEEP);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      ::glDisable(GL_DEPTH_TEST);
      ::glDisable(GL_LIGHTING);
      { // draw shadow on floor
        ::glBegin(GL_QUADS);
        glColor4f(0.1f, 0.1f, 0.1f, 0.5f);
        ::glNormal3d(0,1,0);
        ::glVertex3d(-wfloor,yfloor,-wfloor);
        ::glVertex3d(+wfloor,yfloor,-wfloor);
        ::glVertex3d(+wfloor,yfloor,+wfloor);
        ::glVertex3d(-wfloor,yfloor,+wfloor);
        ::glEnd();
      }
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
      glDisable(GL_STENCIL_TEST);
    }
  }
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int Ptn035[]={   /* # */
   2,  12,66, 87,66,
   2,  12,33, 87,33,
   2,  37,91, 37, 8,
   2,  62, 8, 62,91,
  -1
};
static int PtnA[]={   /* A */
   3,   0, 0,  50,100, 100,0,
   2,  25,50,  75, 50,
  -1,
};
static int PtnB[]={   /* B */
  10,   0, 0, 0,100, 75,100, 87,91, 87,58, 75,50, 100,33, 100,8, 87,0, 0,0,
   2,  75,50, 0, 50,
  -1
};
static int PtnC[]={   /* C */
   8,  100,83, 75,100, 25,100, 0,83, 0,16, 25,0, 75,0, 100,16,
  -1,
};
static int PtnD[]={   /* D */
   7,  0,100, 75,100, 100,83, 100,16, 75,0, 0,0, 0,100,
  -1,
};
static int PtnE[]={   /* E */
   4,  100,100,  0,100, 0,0, 100,0,
   2,    0, 50, 87, 50,
  -1,
};
static int PtnF[]={   /* F */
   3,  100,100,  0,100, 0,0,
   2,    0, 50, 75, 50,
  -1,
};
static int PtnG[]={   /* G */
  10,  100,83, 75,100, 25,100, 0,83, 0,16, 25,0, 75,0, 100,16, 100,41, 62,41,
  -1,
};
///////////
static int Ptn3[]={   /* 3 */
  11,  12,83, 37,100, 75,100, 100,83, 100,66, 75,50, 100,33, 100,16, 75,0, 25,0, 0,16,
  -1
};
static int Ptn4[]={   /* 4 */
   3,  37,100, 12,25, 87,25,
   2,  62, 75, 62, 0,
  -1
};
static int Ptn5[]={   /* 5 */
  10,  87,100, 12,100, 12,41, 37,58, 62,58, 87,41, 87,16, 62,0, 37,0, 12,16,
  -1
};
static int Ptn6[]={   /* 6 */
  12,  87,83, 62,100, 25,100, 0,83, 0,16, 25,0, 75,0, 100,16, 100,33, 75,50, 25,50, 0,33,
  -1
};
static int Ptn7[]={   /* 7 */
   5,  12,83, 12,100, 87,100, 50,33, 50,0,
  -1
};
static int Ptn8[]={   /* 8 */
   9,  100,83, 75,100, 25,100,  0,83,  0,66,  25,50,  75,50, 100,66, 100,83,
   8,   25,50,  0, 33,  0, 16, 25, 0, 75, 0, 100,16, 100,33,  75,50,
  -1
};
static int Ptn9[]={   /* 9 */
  12,  0,16, 25,0, 75,0, 100,16, 100,83, 75,100, 25,100, 0,83, 0,58, 25,41, 75,41, 100,58,
  -1
};

// x = ax*[x] + bx
// y = ay*[y] + by
void DrawCharacter
(int* pChr,
 double ax, double bx,
 double ay, double by)
{
  assert(pChr!=0);
  int icur = 0;
  for(;;){
    int np = pChr[icur];
    if( np == -1 ) break;
    ::glBegin(GL_LINE_STRIP);
    for(int ip=0;ip<np;ip++){
      int ix0 = pChr[icur+1+ip*2+0];
      int iy0 = pChr[icur+1+ip*2+1];
      double sx = ix0*ax+bx;
      double sy = iy0*ay+by;
      ::glVertex2d(sx,sy);
    }
    ::glEnd();
    icur += np*2+1;
  }
}

// x = ax*[x] + bx
// y = ay*[y] + by
void DrawCharacter
(char ic,
 double ax, double bx,
 double ay, double by)
{
  int* pChr = 0;
  if( ic == '#'){ pChr = Ptn035; }
  if( ic == 'A'){ pChr = PtnA; }
  if( ic == 'B'){ pChr = PtnB; }
  if( ic == 'C'){ pChr = PtnC; }
  if( ic == 'D'){ pChr = PtnD; }
  if( ic == 'E'){ pChr = PtnE; }
  if( ic == 'F'){ pChr = PtnF; }
  if( ic == 'G'){ pChr = PtnG; }
  if( ic == '3'){ pChr = Ptn3; }
  if( ic == '4'){ pChr = Ptn4; }
  if( ic == '5'){ pChr = Ptn5; }
  if( ic == '6'){ pChr = Ptn6; }
  if( ic == '7'){ pChr = Ptn7; }
  if( ic == '8'){ pChr = Ptn8; }
  if( ic == '9'){ pChr = Ptn9; }  
  assert(pChr!=0);
  int icur = 0;
  for(;;){
    int np = pChr[icur];
    if( np == -1 ) break;
    ::glBegin(GL_LINE_STRIP);
    for(int ip=0;ip<np;ip++){
      int ix0 = pChr[icur+1+ip*2+0];
      int iy0 = pChr[icur+1+ip*2+1];
      double sx = ix0*ax+bx;
      double sy = iy0*ay+by;
      ::glVertex2d(sx,sy);
    }
    ::glEnd();
    icur += np*2+1;
  }
}

static void UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3])
{
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const double invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

void myGlVertex3d
(unsigned int ixyz,
 const std::vector<double>& aXYZ )
{
  ::glVertex3d(aXYZ[ixyz*3+0], aXYZ[ixyz*3+1], aXYZ[ixyz*3+2] );
}

void myGlVertex2d
(unsigned int ixy,
 const std::vector<double>& aXY )
{
  ::glVertex2d(aXY[ixy*2+0], aXY[ixy*2+1] );
}


void drawLoop2d
(const std::vector<double>& vec)
{
  ::glBegin(GL_LINES);
  const unsigned int nvec = (int)vec.size()/2;
  for(unsigned int ivec=0;ivec<nvec;ivec++){
    unsigned int jvec = ivec+1; if( jvec >= nvec ){ jvec -= nvec; }
    myGlVertex2d(ivec,vec);
    myGlVertex2d(jvec,vec);
  }
  ::glEnd();
  ////
  ::glBegin(GL_POINTS);
  for(unsigned int ivec=0;ivec<nvec;ivec++){
    myGlVertex2d(ivec,vec);
  }
  ::glEnd();
}


void myGlNorm3d
(unsigned int ixyz,
 const std::vector<double>& aNorm )
{
  ::glNormal3d(aNorm[ixyz*3+0], aNorm[ixyz*3+1], aNorm[ixyz*3+2] );
}

void DrawSingleTri3D_FaceNorm
(const std::vector<double>& aXYZ,
 const int* aIndXYZ,
 const double* pUV)
{
  const int i0 = aIndXYZ[0]; assert( i0>=0&&i0<(int)aXYZ.size()/3 );
  const int i1 = aIndXYZ[1]; assert( i1>=0&&i1<(int)aXYZ.size()/3 );
  const int i2 = aIndXYZ[2]; assert( i2>=0&&i2<(int)aXYZ.size()/3 );
  if( i0 == -1 ){
    assert(i1==-1); assert(i2==-1);
    return;
  }
  const double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
  const double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
  const double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
  double un[3], area; UnitNormalAreaTri3D(un,area, p0,p1,p2);
  ::glNormal3dv(un);
  if( pUV != 0 ){ ::glTexCoord2d(pUV[0],pUV[1]); }
  ::glVertex3dv(p0);
  if( pUV != 0 ){ ::glTexCoord2d(pUV[2],pUV[3]); }
  ::glVertex3dv(p1);
  if( pUV != 0 ){ ::glTexCoord2d(pUV[4],pUV[5]); }
  ::glVertex3dv(p2);
}

void DrawSingleQuad3D_FaceNorm
(const std::vector<double>& aXYZ,
 const int* aIndXYZ,
 const double* pUV)
{
  const int i0 = aIndXYZ[0]; assert( i0 >= 0 && i0 < (int)aXYZ.size()/3 );
  const int i1 = aIndXYZ[1]; assert( i1 >= 0 && i1 < (int)aXYZ.size()/3 );
  const int i2 = aIndXYZ[2]; assert( i2 >= 0 && i2 < (int)aXYZ.size()/3 );
  const int i3 = aIndXYZ[3]; assert( i3 >= 0 && i3 < (int)aXYZ.size()/3 );
  if( i0 == -1 ){
    assert(i1==-1 && i2==-1 && i3 ==-1);
    return;
  }
  double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
  double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
  double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
  double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
  {
    double un0[3], area; UnitNormalAreaTri3D(un0,area,  p0, p1, p3);
    if( pUV != 0 ){ ::glTexCoord2d(pUV[0],pUV[1]); }
    ::glNormal3dv(un0);
    ::glVertex3dv(p0);
  }
  {
    double un1[3], area; UnitNormalAreaTri3D(un1,area,  p0, p1, p2);
    if( pUV != 0 ){ ::glTexCoord2d(pUV[2],pUV[3]); }
    ::glNormal3dv(un1);
    ::glVertex3dv(p1);
  }
  {
    double un2[3], area; UnitNormalAreaTri3D(un2,area,  p1, p2, p3);
    if( pUV != 0 ){ ::glTexCoord2d(pUV[4],pUV[5]); }
    ::glNormal3dv(un2);
    ::glVertex3dv(p2);
  }
  {
    double un3[3], area; UnitNormalAreaTri3D(un3,area,  p2, p3, p0);
    if( pUV != 0 ){ ::glTexCoord2d(pUV[6],pUV[7]); }
    ::glNormal3dv(un3);
    ::glVertex3dv(p3);
  }
}

void DrawMeshTri3D_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri)
{
  const int nTri = (int)aTri.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<nTri;++itri){
    DrawSingleTri3D_FaceNorm(aXYZ, aTri.data()+itri*3,0);
  }
  ::glEnd();
}

void DrawMeshElem3D_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem)
{
  const int nelem = aElemInd.size()-1;
  for(int ielem=0;ielem<nelem;++ielem){
    const int ielemind0 = aElemInd[ielem];
    const int ielemind1 = aElemInd[ielem+1];
    if( ielemind1 - ielemind0 == 3 ){
      ::glBegin(GL_TRIANGLES);
      DrawSingleTri3D_FaceNorm(aXYZ, aElem.data()+ielemind0,0);
      ::glEnd();
    }
    else if(ielemind1-ielemind0 == 4){
      ::glBegin(GL_QUADS);
      DrawSingleQuad3D_FaceNorm(aXYZ,aElem.data()+ielemind0,0);
      ::glEnd();
    }
  }
}

void DrawMeshElem3D_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 const std::vector<double>& aUV)
{
  const int nelem = aElemInd.size()-1;
  for(int ielem=0;ielem<nelem;++ielem){
    const int ielemind0 = aElemInd[ielem];
    const int ielemind1 = aElemInd[ielem+1];
    if( ielemind1 - ielemind0 == 3 ){
      ::glBegin(GL_TRIANGLES);
      DrawSingleTri3D_FaceNorm(aXYZ,
                               aElem.data()+ielemind0,
                               aUV.data()+ielemind0*2);
      ::glEnd();
    }
    else if(ielemind1-ielemind0 == 4){
      ::glBegin(GL_QUADS);
      DrawSingleQuad3D_FaceNorm(aXYZ,
                                aElem.data()+ielemind0,
                                aUV.data()+ielemind0*2);
      ::glEnd();
    }
  }
}

void DrawMeshElemPart3D_FaceNorm_TexPoEl
(const std::vector<double>& aXYZ,
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 const std::vector<int>& aIndElemPart,
 const std::vector<double>& aUV)
{
  const bool isUV = (aUV.size()==aElem.size()*2);
  for(int iie=0;iie<(int)aIndElemPart.size();++iie){
    const int ielem = aIndElemPart[iie];
    const int ielemind0 = aElemInd[ielem];
    const int ielemind1 = aElemInd[ielem+1];
    const double* pUV = isUV ? aUV.data()+ielemind0*2:0;
    if( ielemind1 - ielemind0 == 3 ){
      ::glBegin(GL_TRIANGLES);
      DrawSingleTri3D_FaceNorm(aXYZ,
                               aElem.data()+ielemind0,
                               pUV);
      ::glEnd();
    }
    else if(ielemind1-ielemind0 == 4){
      ::glBegin(GL_QUADS);
      DrawSingleQuad3D_FaceNorm(aXYZ,
                                aElem.data()+ielemind0,
                                pUV);
      ::glEnd();
    }
  }
}

void DrawMeshTri3D_FaceNorm_XYsym
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri)
{
  const int nTri = (int)aTri.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<nTri;++itri){
    const int i1 = aTri[itri*3+0];
    const int i2 = aTri[itri*3+1];
    const int i3 = aTri[itri*3+2];
    if( i1 == -1 ){
      assert(i2==-1); assert(i3==-1);
      continue;
    }
    assert( i1 >= 0 && i1 < (int)aXYZ.size()/3 );
    assert( i2 >= 0 && i2 < (int)aXYZ.size()/3 );
    assert( i3 >= 0 && i3 < (int)aXYZ.size()/3 );
    double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], -aXYZ[i1*3+2]};
    double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], -aXYZ[i2*3+2]};
    double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], -aXYZ[i3*3+2]};
    double un[3], area;
    UnitNormalAreaTri3D(un,area, p1,p3,p2);
    ::glNormal3dv(un);
    ::glVertex3dv(p1);
    ::glVertex3dv(p3);
    ::glVertex3dv(p2);
  }
  ::glEnd();
}


void DrawMeshTri3DPart_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<int>& aIndTri)
{
  ::glBegin(GL_TRIANGLES);
  for(int iitri=0;iitri<(int)aIndTri.size();++iitri){
    const int itri = aIndTri[iitri];
    assert( itri>=0&&itri<(int)aTri.size()/3 );
    DrawSingleTri3D_FaceNorm(aXYZ, aTri.data()+itri*3,0);
  }
  ::glEnd();
}

void DrawMeshTri3D_FaceNorm_Flg
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 int iflg,
 const std::vector<int>& aFlgTri)
{
  const int nTri = (int)aTri.size()/3;
//  const int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<nTri;++itri){
    const int iflg0 = aFlgTri[itri];
    if( iflg0 != iflg ) continue;
    DrawSingleTri3D_FaceNorm(aXYZ, aTri.data()+itri*3,0);
  }
  ::glEnd();
}

void DrawMeshTri3D_FaceNormEdge
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri)
{
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ////
  const int nTri = (int)aTri.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for (int itri=0; itri<nTri; ++itri){
    DrawSingleTri3D_FaceNorm(aXYZ, aTri.data()+itri*3,0);
  }
  ::glEnd();

  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (int itri = 0; itri<nTri; ++itri){
    const int i1 = aTri[itri*3+0];
    const int i2 = aTri[itri*3+1];
    const int i3 = aTri[itri*3+2];
    if( i1 == -1 ){
      assert(i2==-1); assert(i3==-1);
      continue;
    }
    assert(i1>=0&&i1 < (int)aXYZ.size()/3 );
    assert(i2>=0&&i2 < (int)aXYZ.size()/3 );
    assert(i3>=0&&i3 < (int)aXYZ.size()/3 );
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glVertex3dv(p3); ::glVertex3dv(p1);
  }
  ::glEnd();

  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}


void DrawMeshTri3D_Edge
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri)
{
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ////
  const int nTri = (int)aTri.size()/3;
  const int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (int itri = 0; itri<nTri; ++itri){
    const int i1 = aTri[itri*3+0];
    const int i2 = aTri[itri*3+1];
    const int i3 = aTri[itri*3+2];
    if( i1 == -1 ){
      assert(i2==-1); assert(i3==-1);
      continue;
    }
    assert(i1>=0&&i1 < nXYZ);
    assert(i2>=0&&i2 < nXYZ);
    assert(i3>=0&&i3 < nXYZ);
    double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glVertex3dv(p3); ::glVertex3dv(p1);
  }
  ::glEnd();

  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}

void DrawMeshTriMap3D_Edge
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<int>& map)
{
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ////
  const int nTri = (int)aTri.size()/3;
  const int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (int itri = 0; itri<nTri; ++itri){
    const int j1 = aTri[itri*3+0];
    const int j2 = aTri[itri*3+1];
    const int j3 = aTri[itri*3+2];
    if( j1 == -1 ){
      assert(j2==-1); assert(j3==-1);
      continue;
    }
    const int i1 = map[j1];
    const int i2 = map[j2];
    const int i3 = map[j3];
    assert(i1>=0 && i1<nXYZ);
    assert(i2>=0 && i2<nXYZ);
    assert(i3>=0 && i3<nXYZ);
    double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glVertex3dv(p3); ::glVertex3dv(p1);
  }
  ::glEnd();
  
  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}

void DrawMeshTri3D_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<double>& aNorm)
{
  unsigned int nTri = (int)aTri.size()/3;
  //  unsigned int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int i1 = aTri[itri*3+0];
    const unsigned int i2 = aTri[itri*3+1];
    const unsigned int i3 = aTri[itri*3+2];
    ::myGlNorm3d(i1,aNorm);  ::myGlVertex3d(i1,aXYZ);
    ::myGlNorm3d(i2,aNorm);  ::myGlVertex3d(i2,aXYZ);
    ::myGlNorm3d(i3,aNorm);  ::myGlVertex3d(i3,aXYZ);
  }
  ::glEnd();
}

void DrawMeshTri3D_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTriVtx,
 const std::vector<double>& aNorm,
 const std::vector<int>& aTriNrm)
{
  const int nTri = (int)aTriVtx.size()/3;
  assert( (int)aTriNrm.size() == nTri*3 );
  //  unsigned int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<nTri;itri++){
    const int iv1 = aTriVtx[itri*3+0];
    const int iv2 = aTriVtx[itri*3+1];
    const int iv3 = aTriVtx[itri*3+2];
    const int in1 = aTriNrm[itri*3+0];
    const int in2 = aTriNrm[itri*3+1];
    const int in3 = aTriNrm[itri*3+2];
    ::myGlNorm3d(in1,aNorm);  ::myGlVertex3d(iv1,aXYZ);
    ::myGlNorm3d(in2,aNorm);  ::myGlVertex3d(iv2,aXYZ);
    ::myGlNorm3d(in3,aNorm);  ::myGlVertex3d(iv3,aXYZ);
  }
  ::glEnd();
}

void DrawMeshTri3D_FaceEdge
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTri)
{
  const int nTri = (int)aTri.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<nTri;itri++){
    const int i1 = aTri[itri*3+0];
    const int i2 = aTri[itri*3+1];
    const int i3 = aTri[itri*3+2];
    ::myGlVertex3d(i1,aXYZ);
    ::myGlVertex3d(i2,aXYZ);
    ::myGlVertex3d(i3,aXYZ);
  }
  ::glEnd();
  ////////
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(int itri=0;itri<nTri;itri++){
    const unsigned int i1 = aTri[itri*3+0];
    const unsigned int i2 = aTri[itri*3+1];
    const unsigned int i3 = aTri[itri*3+2];
    ::myGlVertex3d(i1,aXYZ);
    ::myGlVertex3d(i2,aXYZ);
    ::myGlVertex3d(i2,aXYZ);
    ::myGlVertex3d(i3,aXYZ);
    ::myGlVertex3d(i3,aXYZ);
    ::myGlVertex3d(i1,aXYZ);
  }
  ::glEnd();
}

/////////////////////////////////////

void DrawMeshTri2D_FaceDisp2D
(std::vector<int>& aTri,
 std::vector<double>& aXY,
 std::vector<double>& aDisp)
{
  const int ntri = (int)aTri.size()/3;
  //  const int nxys = (int)aXY.size()/2;
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(int itri=0;itri<ntri;itri++){
    //      double color[3];
    const int i0 = aTri[itri*3+0];
    const int i1 = aTri[itri*3+1];
    const int i2 = aTri[itri*3+2];
    const double p0[2] = { aXY[i0*2+0]+aDisp[i0*2+0], aXY[i0*2+1]+aDisp[i0*2+1] };
    const double p1[2] = { aXY[i1*2+0]+aDisp[i1*2+0], aXY[i1*2+1]+aDisp[i1*2+1] };
    const double p2[2] = { aXY[i2*2+0]+aDisp[i2*2+0], aXY[i2*2+1]+aDisp[i2*2+1] };
    ::glVertex2dv( p0 );
    ::glVertex2dv( p1 );
    ::glVertex2dv( p2 );
  }
  ::glEnd();
  ////////////////
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(int itri=0;itri<ntri;itri++){
    const int i0 = aTri[itri*3+0];
    const int i1 = aTri[itri*3+1];
    const int i2 = aTri[itri*3+2];
    const double p0[2] = { aXY[i0*2+0]+aDisp[i0*2+0], aXY[i0*2+1]+aDisp[i0*2+1] };
    const double p1[2] = { aXY[i1*2+0]+aDisp[i1*2+0], aXY[i1*2+1]+aDisp[i1*2+1] };
    const double p2[2] = { aXY[i2*2+0]+aDisp[i2*2+0], aXY[i2*2+1]+aDisp[i2*2+1] };
    ::glVertex2dv( p0 ); ::glVertex2dv( p1 );
    ::glVertex2dv( p1 ); ::glVertex2dv( p2 );
    ::glVertex2dv( p2 ); ::glVertex2dv( p0 );
  }
  ::glEnd();
}


void DrawMeshTri2D_Edge
(std::vector<int>& aTri,
 std::vector<double>& aXY)
{
  const unsigned int ntri = (int)aTri.size()/3;
  //  const unsigned int nxys = (int)aXY.size()/2;
  ////////////////
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itri=0;itri<ntri;itri++){
    const unsigned int ino0 = aTri[itri*3+0];
    const unsigned int ino1 = aTri[itri*3+1];
    const unsigned int ino2 = aTri[itri*3+2];
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
  }
  ::glEnd();
}

void DrawPoints2D_Vectors
(std::vector<double>& aXY,
 std::vector<double>& aVal,
 int nstride,
 int noffset,
 double mag)
{
  //  const int ntri = (int)aTri.size()/3;
  const int nxys = (int)aXY.size()/2;
  ////
  ::glBegin(GL_LINES);
  for(int ino=0;ino<nxys;ino++){
    const double vx = aVal[ino*nstride+noffset+0]*mag;
    const double vy = aVal[ino*nstride+noffset+1]*mag;
    const double p0[2] = { aXY[ino*2+0],    aXY[ino*2+1]    };
    const double p1[2] = { aXY[ino*2+0]+vx, aXY[ino*2+1]+vy };
    ::glVertex2dv( p0 );
    ::glVertex2dv( p1 );
  }
  ::glEnd();
}

void DrawPoints2D_Points(std::vector<double>& aXY)
{
  const int nxys = (int)aXY.size()/2;
  ::glBegin(GL_POINTS);
  for(int ino=0;ino<nxys;ino++){
    const double p0[2] = { aXY[ino*2+0], aXY[ino*2+1] };
    ::glVertex2dv( p0 );
  }
  ::glEnd();
}

void DrawPoints3D_Points(std::vector<double>& aXYZ)
{
  const int nxyz = (int)aXYZ.size()/3;
  ::glBegin(GL_POINTS);
  for(int ino=0;ino<nxyz;ino++){
    const double p0[3] = { aXYZ[ino*3+0], aXYZ[ino*3+1], aXYZ[ino*3+2]};
    ::glVertex3dv( p0 );
  }
  ::glEnd();
}

////////////////////////////////////////////////////////////////////////////////

void DrawMeshQuad3D_Edge
(const std::vector<double>& aXYZ,
 const std::vector<int>& aQuad)
{
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ////
  const int nQuad = (int)aQuad.size()/4;
  const int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (int iq = 0; iq<nQuad; ++iq){
    const int i1 = aQuad[iq*4+0];
    const int i2 = aQuad[iq*4+1];
    const int i3 = aQuad[iq*4+2];
    const int i4 = aQuad[iq*4+3];
    if( i1 == -1 ){
      assert(i2==-1); assert(i3==-1);
      continue;
    }
    assert(i1>=0&&i1<nXYZ);
    assert(i2>=0&&i2<nXYZ);
    assert(i3>=0&&i3<nXYZ);
    assert(i4>=0&&i4<nXYZ);
    double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double p4[3] = { aXYZ[i4*3+0], aXYZ[i4*3+1], aXYZ[i4*3+2] };
    ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glVertex3dv(p3); ::glVertex3dv(p4);
    ::glVertex3dv(p4); ::glVertex3dv(p1);
  }
  ::glEnd();
  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}

///////////////////////////////////////////

void DrawMeshTet3DSurface_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTet,
 const std::vector<int>& aTetFace)
{
  const unsigned int noelTetFace[4][3] = {
    { 1, 2, 3 },
    { 0, 3, 2 },
    { 0, 1, 3 },
    { 0, 2, 1 } };
  //  const int nTri = (int)aTri.size()/3;
  // const int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glBegin(GL_TRIANGLES);
  for(int itf=0;itf<(int)aTetFace.size()/2;++itf){
    int itet = aTetFace[itf*2+0];
    int iface = aTetFace[itf*2+1];
    const int i1 = aTet[itet*4+noelTetFace[iface][0]];
    const int i2 = aTet[itet*4+noelTetFace[iface][1]];
    const int i3 = aTet[itet*4+noelTetFace[iface][2]];
    if( i1 == -1 ){
      assert(i2==-1); assert(i3==-1);
      continue;
    }
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
  }
  ::glEnd();
}

void DrawMeshTet3DSurface_Edge
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTet,
 const std::vector<int>& aTetFace)
{
  const unsigned int noelTetFace[4][3] = {
    { 1, 2, 3 },
    { 0, 3, 2 },
    { 0, 1, 3 },
    { 0, 2, 1 } };
  
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  ////
  const int nXYZ = (int)aXYZ.size()/3;
  /////
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (int itf=0; itf<(int)aTetFace.size()/2;++itf){
    int itet = aTetFace[itf*2+0];
    int iface = aTetFace[itf*2+1];
    const int i1 = aTet[itet*4+noelTetFace[iface][0]];
    const int i2 = aTet[itet*4+noelTetFace[iface][1]];
    const int i3 = aTet[itet*4+noelTetFace[iface][2]];
    if( i1 == -1 ){
      assert(i2==-1); assert(i3==-1);
      continue;
    }
    assert(i1>=0&&i1 < nXYZ);
    assert(i2>=0&&i2 < nXYZ);
    assert(i3>=0&&i3 < nXYZ);
    double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glVertex3dv(p3); ::glVertex3dv(p1);
  }
  ::glEnd();
  
  if (is_lighting){ ::glEnable(GL_LIGHTING); }
}

void DrawMeshTet3D_Edge
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTet)
{
  for (int itet = 0; itet<(int)aTet.size()/4; itet++){
    const int i0 = aTet[itet*4+0];
    const int i1 = aTet[itet*4+1];
    const int i2 = aTet[itet*4+2];
    const int i3 = aTet[itet*4+3];
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    //::glColor3d(0, 0, 0);
    ::glBegin(GL_LINES);
    ::glVertex3dv(p0); ::glVertex3dv(p1);
    ::glVertex3dv(p0); ::glVertex3dv(p2);
    ::glVertex3dv(p0); ::glVertex3dv(p3);
    ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glVertex3dv(p1); ::glVertex3dv(p3);
    ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glEnd();
  }
}

void DrawMeshTet3D_EdgeDisp
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTet,
 const std::vector<double>& aDisp)
{
  for (int itet = 0; itet<(int)aTet.size()/4; itet++){
    const int i0 = aTet[itet*4+0];
    const int i1 = aTet[itet*4+1];
    const int i2 = aTet[itet*4+2];
    const int i3 = aTet[itet*4+3];
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double q0[3] = { p0[0]+aDisp[i0*3+0], p0[1]+aDisp[i0*3+1], p0[2]+aDisp[i0*3+2] };
    double q1[3] = { p1[0]+aDisp[i1*3+0], p1[1]+aDisp[i1*3+1], p1[2]+aDisp[i1*3+2] };
    double q2[3] = { p2[0]+aDisp[i2*3+0], p2[1]+aDisp[i2*3+1], p2[2]+aDisp[i2*3+2] };
    double q3[3] = { p3[0]+aDisp[i3*3+0], p3[1]+aDisp[i3*3+1], p3[2]+aDisp[i3*3+2] };
    ::glBegin(GL_LINES);
    ::glVertex3dv(q0); ::glVertex3dv(q1);
    ::glVertex3dv(q0); ::glVertex3dv(q2);
    ::glVertex3dv(q0); ::glVertex3dv(q3);
    ::glVertex3dv(q1); ::glVertex3dv(q2);
    ::glVertex3dv(q1); ::glVertex3dv(q3);
    ::glVertex3dv(q2); ::glVertex3dv(q3);
    ::glEnd();
  }
}

void DrawMeshTet3D_FaceNormal
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTet)
{
  for (int itet = 0; itet<(int)aTet.size()/4; itet++){
    const int i0 = aTet[itet*4+0];
    const int i1 = aTet[itet*4+1];
    const int i2 = aTet[itet*4+2];
    const int i3 = aTet[itet*4+3];
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double un0[3], a0; UnitNormalAreaTri3D(un0,a0, p1,p2,p3);
    double un1[3], a1; UnitNormalAreaTri3D(un1,a1, p2,p0,p3);
    double un2[3], a2; UnitNormalAreaTri3D(un2,a2, p3,p0,p1);
    double un3[3], a3; UnitNormalAreaTri3D(un3,a3, p0,p2,p1);
    //    ::glColor3d(0, 0, 0);
    ::glBegin(GL_TRIANGLES);
    ::glNormal3dv(un0); ::glVertex3dv(p1); ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glNormal3dv(un1); ::glVertex3dv(p2); ::glVertex3dv(p3); ::glVertex3dv(p0);
    ::glNormal3dv(un2); ::glVertex3dv(p3); ::glVertex3dv(p0); ::glVertex3dv(p1);
    ::glNormal3dv(un3); ::glVertex3dv(p0); ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glEnd();
  }
}

void DrawHex3D_FaceNorm
(const std::vector<double>& aXYZ,
 const std::vector<int>& aHex)
{
  const int noelElemFace_Hex[8][4] = { // this is corresponds to VTK_VOXEL
    { 0, 4, 7, 3 }, // -x
    { 1, 2, 6, 5 }, // +x
    { 0, 1, 5, 4 }, // -y
    { 3, 7, 6, 2 }, // +y
    { 0, 3, 2, 1 }, // -z
    { 4, 5, 6, 7 }  // +z
  };
  ::glBegin(GL_TRIANGLES);
  for (int ihex = 0; ihex<(int)aHex.size()/8; ihex++){
    const int i0 = aHex[ihex*8+0];
    const int i1 = aHex[ihex*8+1];
    const int i2 = aHex[ihex*8+2];
    const int i3 = aHex[ihex*8+3];
    const int i4 = aHex[ihex*8+4];
    const int i5 = aHex[ihex*8+5];
    const int i6 = aHex[ihex*8+6];
    const int i7 = aHex[ihex*8+7];
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    const double p4[3] = { aXYZ[i4*3+0], aXYZ[i4*3+1], aXYZ[i4*3+2] };
    const double p5[3] = { aXYZ[i5*3+0], aXYZ[i5*3+1], aXYZ[i5*3+2] };
    const double p6[3] = { aXYZ[i6*3+0], aXYZ[i6*3+1], aXYZ[i6*3+2] };
    const double p7[3] = { aXYZ[i7*3+0], aXYZ[i7*3+1], aXYZ[i7*3+2] };
    const double* aP[8] = {p0,p1,p2,p3,p4,p5,p6,p7};
    for(int iface=0;iface<6;++iface){
      const double* q0 = aP[ noelElemFace_Hex[iface][0] ];
      const double* q1 = aP[ noelElemFace_Hex[iface][1] ];
      const double* q2 = aP[ noelElemFace_Hex[iface][2] ];
      const double* q3 = aP[ noelElemFace_Hex[iface][3] ];
      double un0[3], a0; UnitNormalAreaTri3D(un0,a0, q0,q1,q2);
      ::glNormal3dv(un0); ::glVertex3dv(q0); ::glVertex3dv(q1); ::glVertex3dv(q2);
      double un1[3], a1; UnitNormalAreaTri3D(un1,a1, q0,q2,q3);
      ::glNormal3dv(un1); ::glVertex3dv(q0); ::glVertex3dv(q2); ::glVertex3dv(q3);
    }
  }
  ::glEnd();
}

void DrawHex3D_FaceNormDirp
(const std::vector<double>& aXYZ,
 const std::vector<int>& aHex,
 const std::vector<double>& aDisp)
{
  const int noelElemFace_Hex[8][4] = { // this is corresponds to VTK_VOXEL
    { 0, 4, 7, 3 }, // -x
    { 1, 2, 6, 5 }, // +x
    { 0, 1, 5, 4 }, // -y
    { 3, 7, 6, 2 }, // +y
    { 0, 3, 2, 1 }, // -z
    { 4, 5, 6, 7 }  // +z
  };
  ::glBegin(GL_TRIANGLES);
  for (int ihex = 0; ihex<(int)aHex.size()/8; ihex++){
    const int i0 = aHex[ihex*8+0];
    const int i1 = aHex[ihex*8+1];
    const int i2 = aHex[ihex*8+2];
    const int i3 = aHex[ihex*8+3];
    const int i4 = aHex[ihex*8+4];
    const int i5 = aHex[ihex*8+5];
    const int i6 = aHex[ihex*8+6];
    const int i7 = aHex[ihex*8+7];
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    const double p4[3] = { aXYZ[i4*3+0], aXYZ[i4*3+1], aXYZ[i4*3+2] };
    const double p5[3] = { aXYZ[i5*3+0], aXYZ[i5*3+1], aXYZ[i5*3+2] };
    const double p6[3] = { aXYZ[i6*3+0], aXYZ[i6*3+1], aXYZ[i6*3+2] };
    const double p7[3] = { aXYZ[i7*3+0], aXYZ[i7*3+1], aXYZ[i7*3+2] };
    ////
    const double r0[3] = { p0[0]+aDisp[i0*3+0], p0[1]+aDisp[i0*3+1], p0[2]+aDisp[i0*3+2] };
    const double r1[3] = { p1[0]+aDisp[i1*3+0], p1[1]+aDisp[i1*3+1], p1[2]+aDisp[i1*3+2] };
    const double r2[3] = { p2[0]+aDisp[i2*3+0], p2[1]+aDisp[i2*3+1], p2[2]+aDisp[i2*3+2] };
    const double r3[3] = { p3[0]+aDisp[i3*3+0], p3[1]+aDisp[i3*3+1], p3[2]+aDisp[i3*3+2] };
    const double r4[3] = { p4[0]+aDisp[i4*3+0], p4[1]+aDisp[i4*3+1], p4[2]+aDisp[i4*3+2] };
    const double r5[3] = { p5[0]+aDisp[i5*3+0], p5[1]+aDisp[i5*3+1], p5[2]+aDisp[i5*3+2] };
    const double r6[3] = { p6[0]+aDisp[i6*3+0], p6[1]+aDisp[i6*3+1], p6[2]+aDisp[i6*3+2] };
    const double r7[3] = { p7[0]+aDisp[i7*3+0], p7[1]+aDisp[i7*3+1], p7[2]+aDisp[i7*3+2] };
    const double* aR[8] = {r0,r1,r2,r3,r4,r5,r6,r7};
    for(int iface=0;iface<6;++iface){
      const double* q0 = aR[ noelElemFace_Hex[iface][0] ];
      const double* q1 = aR[ noelElemFace_Hex[iface][1] ];
      const double* q2 = aR[ noelElemFace_Hex[iface][2] ];
      const double* q3 = aR[ noelElemFace_Hex[iface][3] ];
      double un0[3], a0; UnitNormalAreaTri3D(un0,a0, q0,q1,q2);
      ::glNormal3dv(un0); ::glVertex3dv(q0); ::glVertex3dv(q1); ::glVertex3dv(q2);
      double un1[3], a1; UnitNormalAreaTri3D(un1,a1, q0,q2,q3);
      ::glNormal3dv(un1); ::glVertex3dv(q0); ::glVertex3dv(q2); ::glVertex3dv(q3);
    }
  }
  ::glEnd();
}

void DrawHex3D_Edge
(const std::vector<double>& aXYZ,
 const std::vector<int>& aHex)
{
  const int noelEdge_Hex[12][2] = {
    {0,1},{3,2},{4,5},{7,6},
    {0,3},{1,2},{4,7},{5,6},
    {0,4},{1,5},{3,7},{2,6} };
  ::glBegin(GL_LINES);
  for (int ihex = 0; ihex<(int)aHex.size()/8; ihex++){
    const int i0 = aHex[ihex*8+0];
    const int i1 = aHex[ihex*8+1];
    const int i2 = aHex[ihex*8+2];
    const int i3 = aHex[ihex*8+3];
    const int i4 = aHex[ihex*8+4];
    const int i5 = aHex[ihex*8+5];
    const int i6 = aHex[ihex*8+6];
    const int i7 = aHex[ihex*8+7];
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    const double p4[3] = { aXYZ[i4*3+0], aXYZ[i4*3+1], aXYZ[i4*3+2] };
    const double p5[3] = { aXYZ[i5*3+0], aXYZ[i5*3+1], aXYZ[i5*3+2] };
    const double p6[3] = { aXYZ[i6*3+0], aXYZ[i6*3+1], aXYZ[i6*3+2] };
    const double p7[3] = { aXYZ[i7*3+0], aXYZ[i7*3+1], aXYZ[i7*3+2] };
    const double* aP[8] = {p0,p1,p2,p3,p4,p5,p6,p7};
    for(int iedge=0;iedge<12;++iedge){
      const double* q0 = aP[ noelEdge_Hex[iedge][0] ];
      const double* q1 = aP[ noelEdge_Hex[iedge][1] ];
      ::glVertex3dv(q0);
      ::glVertex3dv(q1);
    }
  }
  ::glEnd();
}

void DrawMeshTet3D_FaceNormDisp
(const std::vector<double>& aXYZ,
 const std::vector<int>& aTet,
 const std::vector<double>& aDisp)
{
  for (int itet = 0; itet<(int)aTet.size()/4; itet++){
    const int i0 = aTet[itet*4+0];
    const int i1 = aTet[itet*4+1];
    const int i2 = aTet[itet*4+2];
    const int i3 = aTet[itet*4+3];
    const double p0[3] = { aXYZ[i0*3+0]+aDisp[i0*3+0], aXYZ[i0*3+1]+aDisp[i0*3+1], aXYZ[i0*3+2]+aDisp[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0]+aDisp[i1*3+0], aXYZ[i1*3+1]+aDisp[i1*3+1], aXYZ[i1*3+2]+aDisp[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0]+aDisp[i2*3+0], aXYZ[i2*3+1]+aDisp[i2*3+1], aXYZ[i2*3+2]+aDisp[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0]+aDisp[i3*3+0], aXYZ[i3*3+1]+aDisp[i3*3+1], aXYZ[i3*3+2]+aDisp[i3*3+2] };
    double un0[3], a0; UnitNormalAreaTri3D(un0,a0, p1,p2,p3);
    double un1[3], a1; UnitNormalAreaTri3D(un1,a1, p2,p0,p3);
    double un2[3], a2; UnitNormalAreaTri3D(un2,a2, p3,p0,p1);
    double un3[3], a3; UnitNormalAreaTri3D(un3,a3, p0,p2,p1);
    //    ::glColor3d(0, 0, 0);
    ::glBegin(GL_TRIANGLES);
    ::glNormal3dv(un0); ::glVertex3dv(p1); ::glVertex3dv(p2); ::glVertex3dv(p3);
    ::glNormal3dv(un1); ::glVertex3dv(p2); ::glVertex3dv(p3); ::glVertex3dv(p0);
    ::glNormal3dv(un2); ::glVertex3dv(p3); ::glVertex3dv(p0); ::glVertex3dv(p1);
    ::glNormal3dv(un3); ::glVertex3dv(p0); ::glVertex3dv(p1); ::glVertex3dv(p2);
    ::glEnd();
  }
}








///////////////////////////////////



void Draw_SurfaceMeshEdge
(unsigned int nXYZ, const double* paXYZ,
 unsigned int nTri, const unsigned int* paTri)
{
  ::glEnableClientState(GL_VERTEX_ARRAY);
  ::glVertexPointer(3 , GL_DOUBLE , 0 , paXYZ);
  ::glBegin(GL_LINES);
  for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int i1 = paTri[itri*3+0];
    const unsigned int i2 = paTri[itri*3+1];
    const unsigned int i3 = paTri[itri*3+2];
    glArrayElement(i1);
    glArrayElement(i2);
    glArrayElement(i2);
    glArrayElement(i3);
    glArrayElement(i3);
    glArrayElement(i1);
  }
  ::glEnd();
  ::glDisableClientState(GL_VERTEX_ARRAY);
  return;
}

void Draw_SurfaceMeshFace
(unsigned int nXYZ, const double* paXYZ,
 unsigned int nTri, const unsigned int* paTri)
{
  ::glEnableClientState(GL_VERTEX_ARRAY);
  ::glVertexPointer(3 , GL_DOUBLE , 0 , paXYZ);
  ::glDrawElements(GL_TRIANGLES , nTri*3 , GL_UNSIGNED_INT , paTri);
  ::glDisableClientState(GL_VERTEX_ARRAY);
  return;
  /*
  /////
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int i1 = paTri[itri*3+0];
    const unsigned int i2 = paTri[itri*3+1];
    const unsigned int i3 = paTri[itri*3+2];
    ::glVertex3dv(paXYZ+i1*3);
    ::glVertex3dv(paXYZ+i2*3);
    ::glVertex3dv(paXYZ+i3*3);
  }
  ::glEnd();
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int i1 = paTri[itri*3+0];
    const unsigned int i2 = paTri[itri*3+1];
    const unsigned int i3 = paTri[itri*3+2];
    ::glVertex3dv(paXYZ+i1*3);
    ::glVertex3dv(paXYZ+i2*3);
    ::glVertex3dv(paXYZ+i2*3);
    ::glVertex3dv(paXYZ+i3*3);
    ::glVertex3dv(paXYZ+i3*3);
    ::glVertex3dv(paXYZ+i1*3);
  }
  ::glEnd();
  */
}


void DrawBox_MinMaxXYZ
(double x_min, double x_max,
 double y_min, double y_max,
 double z_min, double z_max)
{// show bounding box
  ::glBegin(GL_LINES);
  ::glVertex3d(x_max,y_min,z_max); ::glVertex3d(x_min,y_min,z_max);
  ::glVertex3d(x_max,y_min,z_min); ::glVertex3d(x_min,y_min,z_min);
  ::glVertex3d(x_max,y_max,z_max); ::glVertex3d(x_min,y_max,z_max);
  ::glVertex3d(x_max,y_max,z_min); ::glVertex3d(x_min,y_max,z_min);
  ::glVertex3d(x_max,y_min,z_max); ::glVertex3d(x_max,y_max,z_max);
  ::glVertex3d(x_min,y_min,z_max); ::glVertex3d(x_min,y_max,z_max);
  ::glVertex3d(x_max,y_min,z_min); ::glVertex3d(x_max,y_max,z_min);
  ::glVertex3d(x_min,y_min,z_min); ::glVertex3d(x_min,y_max,z_min);
  ::glVertex3d(x_max,y_min,z_min); ::glVertex3d(x_max,y_min,z_max);
  ::glVertex3d(x_min,y_min,z_min); ::glVertex3d(x_min,y_min,z_max);
  ::glVertex3d(x_max,y_max,z_min); ::glVertex3d(x_max,y_max,z_max);
  ::glVertex3d(x_min,y_max,z_min); ::glVertex3d(x_min,y_max,z_max);
  ::glEnd();
}

void DrawBox_MinMaxXYZ
(double aabbMinMaxXYZ[6])
{// show bounding box
  DrawBox_MinMaxXYZ(aabbMinMaxXYZ[0], aabbMinMaxXYZ[1],
                    aabbMinMaxXYZ[2], aabbMinMaxXYZ[3],
                    aabbMinMaxXYZ[4], aabbMinMaxXYZ[5]);
}


void SaveImage(const std::string& path)
{
  static unsigned int inum = 0;
  int viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  void* image = malloc(3*viewport[2]*viewport[3]);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, viewport[2], viewport[3], GL_RGB, GL_UNSIGNED_BYTE, image);
  unsigned int width = viewport[2];
  unsigned int height = viewport[3];
  std::ofstream fout;
  //  //  fname << "out";
  fout.open(path.c_str(), std::ios::out);
  fout<<"P3\n";
  fout<<width<<" "<<height<<"\n";
  fout<<"255\n";
  //  fout << "255\n";
  //  fout << "255\n";
  char* img = (char*)image;
  for (unsigned int ih = 0; ih<height; ih++){
    for (unsigned int iw = 0; iw<width; iw++){
      unsigned int i = (height-1-ih)*width+iw;
      int r = (unsigned char)img[i*3+0];
      int g = (unsigned char)img[i*3+1];
      int b = (unsigned char)img[i*3+2];
      fout<<r<<" "<<g<<" "<<b<<"\n";
      //    std::cout << i << " " << r << " "<< g << " "<< b << std::endl;
    }
  }
  fout.close();
  //  if( inum >= 700 ) abort();
  //  if( inum >= 400 ) abort();
  if (inum>=600) abort();
  inum++;
}


void LoadImage_PPM
(const std::string& filename,
 std::vector<unsigned char>& image,
 int& width, int& height)
{
  std::ifstream file(filename.c_str(), std::ios::binary);
  if (!file) {
    std::cerr<<"Could not open file \""<<filename<<"\"."<<std::endl;
    return;
  }
  
  std::string header;
  {
    char buff[256];
    file.getline(buff, 256);
    header = std::string(buff);
  }
  if (header=="P6") {
    {
      int max;
      char buff[256];
      file.getline(buff,256);
      if (buff[0]=='#'){
        file>>width>>height>>max;
      }
      else{
        std::stringstream ss(buff);
        ss>>width>>height>>max;
      }
      //      std::cout<<header<<" "<<width<<" "<<height<<" "<<max<<std::endl;
    }
    // max is supporse to be 255
    const int size = width*height*3+1;
    std::vector<unsigned char> data(size);
    file.read(reinterpret_cast<char*>(&(data[0])), size);
    image.resize(3*height*width+256);
    for (int row = 0; row < height; ++row) {
      for (int col = 0; col < width; ++col) {
        int dest_index = 3*((height-row-1) * width+col);
        int src_index = (row * width+col)*3;
        image[dest_index+0] = data[src_index+1];
        image[dest_index+1] = data[src_index+2];
        image[dest_index+2] = data[src_index+3];
      }
    }
  }
  file.close();
}

int ReadPPM_SetTexture(const std::string& fname)
{
  std::cout << "ReadPPM " << std::endl;
  FILE* fp = fopen(fname.c_str(),"r");
  if( fp == NULL ){
    std::cout << "Read PPM Fail" << std::endl;
    return -1;
  }
  
  int w, h;
  std::vector<char> aRGB;
  {
    const unsigned int buffSize = 256;
    char buff[buffSize];
    fgets(buff,buffSize,fp); std::cout << buff << std::endl;
    fgets(buff,buffSize,fp); std::cout << buff << std::endl;
    sscanf(buff,"%d%d",&w,&h);
    fgets(buff,buffSize,fp);  // read 255
  }
  std::cout << "tex size : " << w << " " << h << std::endl;
  //  assert( w >= 0 && h >=0 );
  aRGB.resize(w*h*3);
  const unsigned int buffSize = (unsigned int)(4*3*w*1.2);	// ÇøÇÂÇ¡Ç∆ó]ï™ñ⁄Ç…Ç∆Ç¡ÇƒÇ®Ç≠
  char* buff = new char [buffSize];
  int icnt = 0;
  while (icnt<w*h*3) {
    fgets(buff,buffSize,fp);
    char* pCur = buff;
    char* pNxt;
    for(;;){
      //      if(      pCur[0] == ' ' ){ assert(0); }
      if(      pCur[1] == ' ' || pCur[1] == '\n' ){ pCur[1]='\0'; pNxt=pCur+2; }
      else if( pCur[2] == ' ' || pCur[2] == '\n' ){ pCur[2]='\0'; pNxt=pCur+3; }
      else if( pCur[3] == ' ' || pCur[3] == '\n' ){ pCur[3]='\0'; pNxt=pCur+4; }
      //      else{ assert(0); }
      unsigned int val = atoi(pCur);
      unsigned int ih = icnt/(w*3);
      unsigned int iw = icnt-ih*w*3;
      aRGB[(h-ih-1)*w*3+iw] = val;
      icnt++;
      if( pNxt[0] == '\n' || pNxt[0] == '\0') break;
      pCur = pNxt;
    }
  }
  delete[] buff;
  //	this->SetImage(w,h,aRGB);
  
  std::cout << "width height : " << w << " " << h << std::endl;
  
  ////////////////
  
  GLubyte* inputRGB = new GLubyte [w*h*3];
  for(int i=0;i<w*h*3;i++){ inputRGB[i] = aRGB[i]; }
  
  int m_texWidth = 256;
  int m_texHeight = 256;
  //	std::cout << m_texWidth << " " << m_texHight << std::endl;
  
  GLubyte* scaledRGB;
  if( w == m_texWidth && h == m_texHeight ){
    scaledRGB = inputRGB;
  }
  else{
    scaledRGB = new GLubyte [m_texWidth*m_texHeight*3];
    gluScaleImage( GL_RGB, w, h, GL_UNSIGNED_BYTE, inputRGB,
                  m_texWidth, m_texHeight, GL_UNSIGNED_BYTE, scaledRGB );
    delete [] inputRGB;
  }
  
  glEnable(GL_TEXTURE_2D);
  GLuint m_texName = 0;
  glGenTextures(1 , &m_texName);
  glBindTexture(GL_TEXTURE_2D , m_texName);
  glTexImage2D(GL_TEXTURE_2D , 0 , 3 , m_texWidth, m_texHeight,
               0 , GL_RGB , GL_UNSIGNED_BYTE , scaledRGB );
  delete[] scaledRGB;
  
  //	std::cout << m_texName << std::endl;
  
  return (int)m_texName;
}

bool LoadTGAFile
(const char *filename,
 SFile_TGA *tgaFile)
{
  FILE *filePtr;
  unsigned char ucharBad;
  short int sintBad;
  long imageSize;
  int colorMode;
  unsigned char colorSwap;
  
  // Open the TGA file.
  filePtr = fopen(filename, "rb");
  if (filePtr == NULL)
  {
    return false;
  }
  
  // Read the two first bytes we don't need.
  fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
  fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
  
  // Which type of image gets stored in imageTypeCode.
  fread(&tgaFile->imageTypeCode, sizeof(unsigned char), 1, filePtr);
  
  // For our purposes, the type code should be 2 (uncompressed RGB image)
  // or 3 (uncompressed black-and-white images).
  if (tgaFile->imageTypeCode != 2 && tgaFile->imageTypeCode != 3)
  {
    fclose(filePtr);
    return false;
  }
  
  // Read 13 bytes of data we don't need.
  fread(&sintBad, sizeof(short int), 1, filePtr);
  fread(&sintBad, sizeof(short int), 1, filePtr);
  fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
  fread(&sintBad, sizeof(short int), 1, filePtr);
  fread(&sintBad, sizeof(short int), 1, filePtr);
  
  // Read the image's width and height.
  fread(&tgaFile->imageWidth, sizeof(short int), 1, filePtr);
  fread(&tgaFile->imageHeight, sizeof(short int), 1, filePtr);
  
  // Read the bit depth.
  fread(&tgaFile->bitCount, sizeof(unsigned char), 1, filePtr);
  
  // Read one byte of data we don't need.
  fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
  
  // Color mode -> 3 = BGR, 4 = BGRA.
  colorMode = tgaFile->bitCount / 8;
  imageSize = tgaFile->imageWidth * tgaFile->imageHeight * colorMode;
  
  // Allocate memory for the image data.
  tgaFile->imageData = (unsigned char*)malloc(sizeof(unsigned char)*imageSize);
  
  // Read the image data.
  fread(tgaFile->imageData, sizeof(unsigned char), imageSize, filePtr);
  
  // Change from BGR to RGB so OpenGL can read the image data.
  for (int imageIdx = 0; imageIdx < imageSize; imageIdx += colorMode)
  {
    colorSwap = tgaFile->imageData[imageIdx];
    tgaFile->imageData[imageIdx] = tgaFile->imageData[imageIdx + 2];
    tgaFile->imageData[imageIdx + 2] = colorSwap;
  }
  
  fclose(filePtr);
  return true;
}

GLuint LoadTexture
(const unsigned char* image,
 const int width, const int height, const int bpp)
{
  GLuint id_tex; glGenTextures(1, &id_tex);
  glBindTexture(GL_TEXTURE_2D, id_tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  if( bpp == 3 ){
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, (GLsizei)width, (GLsizei)height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
  }
  if( bpp == 4 ){
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, (GLsizei)width, (GLsizei)height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
  }
  return id_tex;
};

void DrawTextureBackground
(const GLuint tex,
 const int imgWidth,
 const int imgHeight,
 const int winWidth,
 const int winHeight)
{
  double imgAsp = (double)imgWidth/imgHeight;
  double winAsp = (double)winWidth/winHeight;
  /////
  glPushAttrib(GL_TRANSFORM_BIT|GL_CURRENT_BIT|GL_ENABLE_BIT);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  
  if (winAsp>imgAsp){
    double imgWidth2 = imgHeight*winAsp;
    gluOrtho2D(
               -0.5*imgWidth2,+0.5*imgWidth2,
               -0.5*imgHeight,+0.5*imgHeight);
  }
  else{
    double imgHeight2 = (double)imgWidth/winAsp;
    gluOrtho2D(
               -0.5*imgWidth, +0.5*imgWidth,
               -0.5*imgHeight2, +0.5*imgHeight2);
  }
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  
  glBindTexture(GL_TEXTURE_2D, tex);
  
  glColor3f(1.0f, 1.0f, 1.0f);
  glBegin(GL_QUADS);
  glTexCoord2i(0, 0); glVertex2d(-0.5*imgWidth, -0.5*imgHeight);
  glTexCoord2i(1, 0); glVertex2d(+0.5*imgWidth, -0.5*imgHeight);
  glTexCoord2i(1, 1); glVertex2d(+0.5*imgWidth, +0.5*imgHeight);
  glTexCoord2i(0, 1); glVertex2d(-0.5*imgWidth, +0.5*imgHeight);
  glEnd();
  
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  
  glPopAttrib();
}

// use it for GLSL shader drawing
void DrawRectangle_FullCanvas()
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  
  glBegin(GL_QUADS);
  glVertex2d(-1, -1);
  glVertex2d(+1, -1);
  glVertex2d(+1, +1);
  glVertex2d(-1, +1);
  glEnd();
  
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  
  glPopAttrib();
}


void CTexManager::Clear(){
  for(int itex=0;itex<(int)aTexInfo.size();++itex){
    unsigned int id_tex_gl = aTexInfo[itex].id_tex_gl;
    if( glIsTexture(id_tex_gl) ){
      ::glDeleteTextures(1, &id_tex_gl);
    }
  }
  aTexInfo.clear();
}

void CTexManager::BindTexturePath(const std::string& path) const {
  for(int itex=0;itex<(int)aTexInfo.size();++itex){
    if( aTexInfo[itex].full_path != path ) continue;
    glBindTexture(GL_TEXTURE_2D, aTexInfo[itex].id_tex_gl );
    glEnable(GL_TEXTURE_2D);
  }
}

void showdepth()
{
  GLint view[4]; glGetIntegerv(GL_VIEWPORT, view); // get viewport size
  GLubyte* buffer = (GLubyte *)malloc(view[2] * view[3]); // get buffer with view port size
  if (!buffer) { return; }
  
  glFinish(); // wait for finishing display
  
  // read depth buffer
  glReadPixels(view[0], view[1], view[2], view[3],
               GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, buffer);
  // write to color buffer
  glDrawPixels(view[2], view[3], GL_LUMINANCE, GL_UNSIGNED_BYTE, buffer);
  
  free(buffer);
}

