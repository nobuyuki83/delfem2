/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/caddtri_v3.h"
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

// -------------------------------------------------


DFM2_INLINE void delfem2::opengl::DrawMeshDynTri_FaceNorm(
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec3d>& aVec3)
{
  //  ::glPushAttrib(GL_ENABLE_BIT);
  
  ::glEnable(GL_LIGHTING);
  //  ::myGlColorDiffuse(CColor::Orange());
  ::glBegin(GL_TRIANGLES);
  for (const auto & itri : aSTri){
    const unsigned int i0 = itri.v[0];
    const unsigned int i1 = itri.v[1];
    const unsigned int i2 = itri.v[2];
    if( i0 == UINT_MAX ){
      assert( i1 == UINT_MAX );
      assert( i2 == UINT_MAX );
      continue;
    }
    {
      CVec3d n; UnitNormal(n, aVec3[i0], aVec3[i1], aVec3[i2]);
      ::glNormal3d(n.x,n.y,n.z);
    }
    {
      CVec3d p0 = aVec3[i0];
      ::glVertex3d(p0.x,p0.y, p0.z);
    }
    {
      CVec3d p1 = aVec3[i1];
      ::glVertex3d(p1.x,p1.y,p1.z);
    }
    {
      CVec3d p2 = aVec3[i2];
      ::glVertex3d(p2.x,p2.y,p2.z);
    }
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshDynTri_FaceNorm(
    const std::vector<CDynTri>& aSTri,
    const double* aXYZ)
{
  //  ::glPushAttrib(GL_ENABLE_BIT);
  
  ::glEnable(GL_LIGHTING);
  //  ::myGlColorDiffuse(CColor::Orange());
  ::glBegin(GL_TRIANGLES);
  for (const auto & tri : aSTri){
    const unsigned int i0 = tri.v[0];
    const unsigned int i1 = tri.v[1];
    const unsigned int i2 = tri.v[2];
    if( i0 == UINT_MAX ){
      assert( i1 == UINT_MAX );
      assert( i2 == UINT_MAX );
      continue;
    }
    const double* p0 = aXYZ+i0*3;
    const double* p1 = aXYZ+i1*3;
    const double* p2 = aXYZ+i2*3;
    {
      double a, n[3];
      UnitNormalAreaTri3(n, a,
                         p0, p1, p2);
      ::glNormal3dv(n);
    }
    ::glVertex3dv(p0);
    ::glVertex3dv(p1);
    ::glVertex3dv(p2);
  }
  ::glEnd();
}


DFM2_INLINE void delfem2::opengl::DrawMeshDynTri_FaceNormTex(
    const std::vector<CDynTri>& aSTri,
    const double* aXYZ,
    const std::vector<CVec2d>& aVec2)
{
  //  ::glPushAttrib(GL_ENABLE_BIT);

  ::glEnable(GL_LIGHTING);
  //  ::myGlColorDiffuse(CColor::Orange());
  ::glBegin(GL_TRIANGLES);
  for (const auto & tri : aSTri){
    const unsigned int i0 = tri.v[0];
    const unsigned int i1 = tri.v[1];
    const unsigned int i2 = tri.v[2];
    if( i0 == UINT_MAX ){
      assert( i1 == UINT_MAX );
      assert( i2 == UINT_MAX );
      continue;
    }
    const double* p0 = aXYZ+i0*3;
    const double* p1 = aXYZ+i1*3;
    const double* p2 = aXYZ+i2*3;
    {
      double a, n[3];
      UnitNormalAreaTri3(n, a,
                         p0, p1, p2);
      ::glNormal3dv(n);
    }
    ::glTexCoord2d(aVec2[i0].x,aVec2[i0].y);
    ::glVertex3dv(p0);
    ::glTexCoord2d(aVec2[i1].x,aVec2[i1].y);
    ::glVertex3dv(p1);
    ::glTexCoord2d(aVec2[i2].x,aVec2[i2].y);
    ::glVertex3dv(p2);
  }
  ::glEnd();
}


DFM2_INLINE void delfem2::opengl::DrawMeshDynTri_Edge(
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec3d>& aVec3)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (const auto & tri : aSTri){
    const unsigned int i0 = tri.v[0];
    const unsigned int i1 = tri.v[1];
    const unsigned int i2 = tri.v[2];
    if( i0 == UINT_MAX ){
      assert( i1 == UINT_MAX );
      assert( i2 == UINT_MAX );
    }
    const CVec3d& p0 = aVec3[i0];
    const CVec3d& p1 = aVec3[i1];
    const CVec3d& p2 = aVec3[i2];
    glVertex3d(p0.x,p0.y,p0.z);
    glVertex3d(p1.x,p1.y,p1.z);
    glVertex3d(p1.x,p1.y,p1.z);
    glVertex3d(p2.x,p2.y,p2.z);
    glVertex3d(p2.x,p2.y,p2.z);
    glVertex3d(p0.x,p0.y,p0.z);
  }
  ::glEnd();
}


DFM2_INLINE void delfem2::opengl::DrawMeshDynTri3D_Edge(
    const std::vector<double>& aXYZ,
    const std::vector<CDynTri>& aSTri)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (const auto & tri : aSTri){
    const unsigned int i0 = tri.v[0];
    const unsigned int i1 = tri.v[1];
    const unsigned int i2 = tri.v[2];
    if( i0 == UINT_MAX ){
      assert( i1 == UINT_MAX );
      assert( i2 == UINT_MAX );
    }
    glVertex3d(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
    glVertex3d(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
    glVertex3d(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
    glVertex3d(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
    glVertex3d(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
    glVertex3d(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
  }
  ::glEnd();
}

// -------------------------------------------------------------------------

