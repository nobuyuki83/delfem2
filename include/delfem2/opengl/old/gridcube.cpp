/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/gridcube.h"
#include "delfem2/srchbv3aabb.h"
#include "delfem2/vec3.h"

#if defined(_WIN32) // windows
  #define NOMINMAX   // to remove min,max macro
  #include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
  #define GL_SILENCE_DEPRECATION
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif

// ------------------------------------

static void myGlNormal(const delfem2::CVec3d& n){ ::glNormal3d(n.x,n.y,n.z); }
static void myGlVertex(const delfem2::CVec3d& v){ ::glVertex3d(v.x,v.y,v.z); }
static void myGlColorDiffuse(float r, float g, float b, float a){
  ::glColor4d(r, g, b, a );
  float c[4] = {r, g, b, a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

// -------------------------------

const unsigned int noelElemFace_Vox[6][4] = {
  { 0, 4, 6, 2 }, // -x
  { 1, 3, 7, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 2, 6, 7, 3 }, // +y
  { 0, 2, 3, 1 }, // -z
  { 4, 5, 7, 6 } }; // +z

const delfem2::CVec3d normalHexFace[6] = {
  delfem2::CVec3d(-1, 0, 0),
  delfem2::CVec3d(+1, 0, 0),
  delfem2::CVec3d( 0,-1, 0),
  delfem2::CVec3d( 0,+1, 0),
  delfem2::CVec3d( 0, 0,-1),
  delfem2::CVec3d( 0, 0,+1)
};

void delfem2::opengl::Draw_CubeGrid(
    bool is_picked,
    int iface_picked,
    double elen,
    const delfem2::CVec3d& org,
    const delfem2::CCubeGrid& cube)
{
  if( !cube.is_active ) return;
  int ih = cube.ivx;
  int jh = cube.ivy;
  int kh = cube.ivz;
  delfem2::CVec3d aP[8] = {
    org + elen * delfem2::CVec3d(ih+0,jh+0,kh+0),
    org + elen * delfem2::CVec3d(ih+1,jh+0,kh+0),
    org + elen * delfem2::CVec3d(ih+0,jh+1,kh+0),
    org + elen * delfem2::CVec3d(ih+1,jh+1,kh+0),
    org + elen * delfem2::CVec3d(ih+0,jh+0,kh+1),
    org + elen * delfem2::CVec3d(ih+1,jh+0,kh+1),
    org + elen * delfem2::CVec3d(ih+0,jh+1,kh+1),
    org + elen * delfem2::CVec3d(ih+1,jh+1,kh+1) };
  ::glEnable(GL_LIGHTING);
  ::glBegin(GL_QUADS);
  for(int iface=0;iface<6;++iface){
    if( is_picked && iface_picked == iface ){
      ::myGlColorDiffuse(1,1,0,1);
    }
    else{
      ::myGlColorDiffuse(0.8f,0.8f,0.8f,1.f);
    }
    myGlNormal(normalHexFace[iface]);
    myGlVertex(aP[noelElemFace_Vox[iface][0]]);
    myGlVertex(aP[noelElemFace_Vox[iface][1]]);
    myGlVertex(aP[noelElemFace_Vox[iface][2]]);
    myGlVertex(aP[noelElemFace_Vox[iface][3]]);
  }
  ::glEnd();
  // --------------------------
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINE_STRIP);
  myGlVertex(aP[0]);
  myGlVertex(aP[1]);
  myGlVertex(aP[3]);
  myGlVertex(aP[2]);
  myGlVertex(aP[0]);
  myGlVertex(aP[4]);
  myGlVertex(aP[5]);
  myGlVertex(aP[7]);
  myGlVertex(aP[6]);
  myGlVertex(aP[4]);
  glEnd();
  ::glBegin(GL_LINES);
  myGlVertex(aP[1]); myGlVertex(aP[5]);
  myGlVertex(aP[2]); myGlVertex(aP[6]);
  myGlVertex(aP[3]); myGlVertex(aP[7]);
  ::glEnd();
}



void Draw(
    const delfem2::CBV3d_AABB& aabb)
{
  double x_min = aabb.bbmin[0];
  double x_max = aabb.bbmax[0];
  double y_min = aabb.bbmin[1];
  double y_max = aabb.bbmax[1];
  double z_min = aabb.bbmin[2];
  double z_max = aabb.bbmax[2];
  //
  const double pxyz[3] = {x_min,y_min,z_min};
  const double pxyZ[3] = {x_min,y_min,z_max};
  const double pxYz[3] = {x_min,y_max,z_min};
  const double pxYZ[3] = {x_min,y_max,z_max};
  const double pXyz[3] = {x_max,y_min,z_min};
  const double pXyZ[3] = {x_max,y_min,z_max};
  const double pXYz[3] = {x_max,y_max,z_min};
  const double pXYZ[3] = {x_max,y_max,z_max};
  ::glColor3d(0,0,0);
  ::glLineWidth(1);
  ::glBegin(GL_LINES);
  ::glVertex3dv(pxyz); ::glVertex3dv(pxyZ);
  ::glVertex3dv(pxYz); ::glVertex3dv(pxYZ);
  ::glVertex3dv(pXyz); ::glVertex3dv(pXyZ);
  ::glVertex3dv(pXYz); ::glVertex3dv(pXYZ);
  ////
  ::glVertex3dv(pxyz); ::glVertex3dv(pXyz);
  ::glVertex3dv(pxyZ); ::glVertex3dv(pXyZ);
  ::glVertex3dv(pxYz); ::glVertex3dv(pXYz);
  ::glVertex3dv(pxYZ); ::glVertex3dv(pXYZ);
  ////
  ::glVertex3dv(pxyz); ::glVertex3dv(pxYz);
  ::glVertex3dv(pxyZ); ::glVertex3dv(pxYZ);
  ::glVertex3dv(pXyz); ::glVertex3dv(pXYz);
  ::glVertex3dv(pXyZ); ::glVertex3dv(pXYZ);
  ::glEnd();
}


