/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>

#include "delfem2/vec3.h"

#if defined(__APPLE__) && defined(__MACH__)
  #include <OpenGL/gl.h>
#elif defined(_WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else
  #include <GL/gl.h>
#endif

#include "delfem2/opengl/glold_funcs.h"
//
#include "delfem2/opengl/caddtri_v3_glold.h"

namespace dfm2 = delfem2;

// -------------------------------------------------


void dfm2::opengl::DrawMeshDynTri_FaceNorm
(const std::vector<CDynTri>& aSTri,
 const std::vector<CVec3d>& aVec3)
{
  //  ::glPushAttrib(GL_ENABLE_BIT);
  
  ::glEnable(GL_LIGHTING);
  //  ::myGlColorDiffuse(CColor::Orange());
  ::glBegin(GL_TRIANGLES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
      continue;
    }
    {
      CVec3d n; UnitNormal(n, aVec3[i0], aVec3[i1], aVec3[i2]);
      ::glNormal3d(n.x(),n.y(),n.z());
    }
    {
      CVec3d p0 = aVec3[i0];
      ::glVertex3d(p0.x(),p0.y(),p0.z());
    }
    {
      CVec3d p1 = aVec3[i1];
      ::glVertex3d(p1.x(),p1.y(),p1.z());
    }
    {
      CVec3d p2 = aVec3[i2];
      ::glVertex3d(p2.x(),p2.y(),p2.z());
    }
  }
  ::glEnd();
}

void dfm2::opengl::DrawMeshDynTri_FaceNorm
 (const std::vector<CDynTri>& aSTri,
  const double* aXYZ)
{
  //  ::glPushAttrib(GL_ENABLE_BIT);
  
  ::glEnable(GL_LIGHTING);
  //  ::myGlColorDiffuse(CColor::Orange());
  ::glBegin(GL_TRIANGLES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
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


void dfm2::opengl::DrawMeshDynTri_Edge
(const std::vector<CDynTri>& aSTri,
 const std::vector<CVec3d>& aVec3)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
    }
    const CVec3d& p0 = aVec3[i0];
    const CVec3d& p1 = aVec3[i1];
    const CVec3d& p2 = aVec3[i2];
    glVertex3d(p0.x(),p0.y(),p0.z());
    glVertex3d(p1.x(),p1.y(),p1.z());
    glVertex3d(p1.x(),p1.y(),p1.z());
    glVertex3d(p2.x(),p2.y(),p2.z());
    glVertex3d(p2.x(),p2.y(),p2.z());
    glVertex3d(p0.x(),p0.y(),p0.z());
  }
  ::glEnd();
}


void dfm2::opengl::DrawMeshDynTri3D_Edge
(const std::vector<double>& aXYZ,
 const std::vector<CDynTri>& aSTri)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
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
