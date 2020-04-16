/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/vec3.h"

#include "glad/glad.h" // gl3.0+
#if defined(__APPLE__) && defined(__MACH__) // Mac
  #include <OpenGL/gl.h>
#elif defined(_WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else
  #include <GL/gl.h>
#endif
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/render2tex_glold.h"

// --------------------------------------------

void delfem2::opengl::CRender2Tex_DrawOldGL::SetPointColor(double r, double g, double b){
  colorPoint[0] = r;
  colorPoint[1] = g;
  colorPoint[2] = b;
}

void delfem2::opengl::CRender2Tex_DrawOldGL::InitGL() {
  CRender2Tex::InitGL();
  if( aRGBA_8ui.size() == nResX*nResY*4 ){
    ::glBindTexture(GL_TEXTURE_2D, id_tex_color);
      // define size and format of level 0
    ::glTexImage2D(GL_TEXTURE_2D, 0,
                   GL_RGBA, nResX, nResY, 0,
                   GL_RGBA, GL_UNSIGNED_BYTE, aRGBA_8ui.data());
  }
}


void delfem2::opengl::CRender2Tex_DrawOldGL::SetView(){
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  delfem2::opengl::ViewTransformation(x_axis,z_axis,origin);
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::glOrtho(0.0, +lengrid*nResX,
            0.0, +lengrid*nResY,
            0, z_range);
  ::glMatrixMode(GL_MODELVIEW);
}

void delfem2::opengl::CRender2Tex_DrawOldGL::Start()
{
  CRender2Tex::Start();
  this->SetView();
}

void delfem2::opengl::CRender2Tex_DrawOldGL::GetDepth()
{
  CRender2Tex::ExtractFromTexture_Depth(aZ);
  for(auto& z: aZ){
    z *= -this->z_range;
  }
}

void delfem2::opengl::CRender2Tex_DrawOldGL::GetColor()
{
  if( is_rgba_8ui ){
    CRender2Tex::ExtractFromTexture_RGBA8UI(aRGBA_8ui);
  }
  else{
    CRender2Tex::ExtractFromTexture_RGBA32F(aRGBA_32f);
  }
}

void delfem2::opengl::CRender2Tex_DrawOldGL::Draw() const {
  ::glPointSize(this->pointSize);
  this->Draw_Point();
  // -----------
  ::glLineWidth(3);
  this->Draw_Axis();
  // ----------
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  this->Draw_BoundingBox();
  // -----------
  if( id_tex_color > 0 && this->isDrawTex ){
    const CVec3d& dx = x_axis;
    const CVec3d& dy = CVec3d(z_axis)^dx;
    const double lx = lengrid*nResX;
    const double ly = lengrid*nResY;
    CVec3d p0 = origin;
    CVec3d p1 = p0 + lx*dx;
    CVec3d p2 = p0 + lx*dx + ly*dy;
    CVec3d p3 = p0 + ly*dy;
    ::glEnable(GL_TEXTURE_2D);
    ::glDisable(GL_LIGHTING);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_color);
    ::glColor3d(1,1,1);
    ::glBegin(GL_QUADS);
    ::glTexCoord2d(0.0, 0.0); opengl::myGlVertex(p0);
    ::glTexCoord2d(1.0, 0.0); opengl::myGlVertex(p1);
    ::glTexCoord2d(1.0, 1.0); opengl::myGlVertex(p2);
    ::glTexCoord2d(0.0, 1.0); opengl::myGlVertex(p3);
    ::glEnd();
    ::glBindTexture(GL_TEXTURE_2D, 0);
    ::glDisable(GL_TEXTURE_2D);
  }
}

void delfem2::opengl::CRender2Tex_DrawOldGL::Draw_Axis() const
{
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  delfem2::opengl::ModelTransformation(x_axis, z_axis, origin);
  delfem2::opengl::DrawAxis(draw_len_axis);
    ::glPopMatrix();
}

void delfem2::opengl::CRender2Tex_DrawOldGL::Draw_BoundingBox() const
{
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  delfem2::opengl::ModelTransformation(x_axis, z_axis, origin);
  ::glLineWidth(3);
  delfem2::opengl::Draw_AABB3D_MinMaxXYZ_Edge(0.0, lengrid*nResX, 0.0, lengrid*nResY, 0.0, -z_range);
  ::glPopMatrix();
}

void delfem2::opengl::CRender2Tex_DrawOldGL::Draw_Point() const
{
  ::glDisable(GL_LIGHTING);
  if( aZ.size() != nResX*nResY ) return;
  if( colorPoint.size() == 3 ){ ::glColor3dv(colorPoint.data()); }
  if( colorPoint.size() == 4 ){ ::glColor4dv(colorPoint.data()); }
  //
  const CVec3d& dx = x_axis;
  const CVec3d& dz = z_axis;
  const CVec3d& dy = Cross(dz,dx);
  ::glBegin(GL_POINTS);
  for(unsigned int iy=0;iy<nResY;++iy){
    for(unsigned int ix=0;ix<nResX;++ix){
      double lz = aZ[iy*nResX+ix];
      if( this->isDrawOnlyHitPoints && lz <= -z_range*0.9 ){ continue; }
      double lx = (ix+0.5)*lengrid;
      double ly = (iy+0.5)*lengrid;
      CVec3d vp = lx*dx+ly*dy+lz*dz + CVec3d(origin);
      delfem2::opengl::myGlVertex(vp);
    }
  }
  ::glEnd();
}

std::vector<double> delfem2::opengl::CRender2Tex_DrawOldGL::getGPos(int ix, int iy) const
{
  const CVec3d& dx = x_axis;
  const CVec3d& dz = z_axis;
  const CVec3d& dy = Cross(dz,dx);
  double lz = -aZ[iy*nResX+ix]*z_range;
  double lx = (ix+0.5)*lengrid;
  double ly = (iy+0.5)*lengrid;
  CVec3d vp = lx*dx+ly*dy+lz*dz + CVec3d(origin);
    //  std::vector<double> res;
    //  res.push_back(vp.x());
    //  res.push_back(vp.y());
    //  res.push_back(vp.z());
  return vp.stlvec();
}


// --------------------------------------------------------

void delfem2::opengl::CRender2Tex_DrawOldGL_BOX::Initialize(unsigned int nresX,
                                                         unsigned int nresY,
                                                         unsigned int nresZ,
                                                         double elen)
{
  aSampler.resize(6);
  aSampler[0].SetTextureProperty(nresY, nresZ, true); // +x
  aSampler[0].SetCoord(elen, elen*nresX,
                       CVec3d(+0.5*elen*nresX,-0.5*elen*nresY,-0.5*elen*nresZ).stlvec(),
                       CVec3d(+1,  0, 0).stlvec(),
                       CVec3d( 0, +1, 0).stlvec() );
  aSampler[0].SetPointColor(1.0, 0.0, 0.0);
  //
  aSampler[1].SetTextureProperty(nresY, nresZ, true); // -x
  aSampler[1].SetCoord(elen, elen*nresX,
                       CVec3d(-0.5*elen*nresX,-0.5*elen*nresY,+0.5*elen*nresZ).stlvec(),
                       CVec3d(-1,  0, 0).stlvec(),
                       CVec3d( 0, +1, 0).stlvec() );
  aSampler[1].SetPointColor(1.0, 0.5, 0.5);
  //
  aSampler[2].SetTextureProperty(nresX, nresZ, true); // +y
  aSampler[2].SetCoord(elen, elen*nresY,
                       CVec3d(-0.5*elen*nresX,+0.5*elen*nresY,+0.5*elen*nresZ).stlvec(),
                       CVec3d(0,+1,0).stlvec(),
                       CVec3d(1,+0,0).stlvec() );
  aSampler[2].SetPointColor(0.0, 1.0, 0.0);
  //
  aSampler[3].SetTextureProperty(nresX, nresZ, true); // -y
  aSampler[3].SetCoord(elen, elen*nresY,
                       CVec3d(-0.5*elen*nresX,-0.5*elen*nresY,-0.5*elen*nresZ).stlvec(),
                       CVec3d(0,-1,0).stlvec(),
                       CVec3d(1,+0,0).stlvec() );
  aSampler[3].SetPointColor(0.5, 1.0, 0.5);
  //
  aSampler[4].SetTextureProperty(nresX, nresY, true);
  aSampler[4].SetCoord(elen, elen*nresZ,
                       CVec3d(-0.5*elen*nresX,-0.5*elen*nresY,+0.5*elen*nresZ).stlvec(),
                       CVec3d(0,0,+1).stlvec(),
                       CVec3d(1,0,0).stlvec() );
  aSampler[4].SetPointColor(0.0, 0.0, 1.0);
  //
  aSampler[5].SetTextureProperty(nresX, nresY, true);
  aSampler[5].SetCoord(elen, elen*nresZ,
                       CVec3d(-0.5*elen*nresX,+0.5*elen*nresY,-0.5*elen*nresZ).stlvec(),
                       CVec3d(0,0,-1).stlvec(),
                       CVec3d(1,0,0).stlvec() );
  aSampler[5].SetPointColor(0.5, 0.5, 1.0);
  // ------------------------
  for(auto& smplr : aSampler){
    smplr.draw_len_axis = 0.2;
    smplr.isDrawTex = false;
    smplr.isDrawOnlyHitPoints = true;
  }
}
