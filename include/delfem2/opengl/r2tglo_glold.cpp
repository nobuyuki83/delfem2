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
#include "delfem2/opengl/r2tglo_glold.h"

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
  const double pmin[3] = {0.0, 0.0, 0.0};
  const double pmax[3] = {lengrid*nResX, lengrid*nResY, -z_range};
  delfem2::opengl::DrawBox3_Edge(pmin, pmax);
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

void delfem2::opengl::CRender2Tex_DrawOldGL::getGPos
 (double* p,
  int ix, int iy) const
{
  const CVec3d& dx = x_axis;
  const CVec3d& dz = z_axis;
  const CVec3d& dy = Cross(dz,dx);
  double lz = aZ[iy*nResX+ix];
  double lx = (ix+0.5)*lengrid;
  double ly = (iy+0.5)*lengrid;
  CVec3d vp = lx*dx+ly*dy+lz*dz + CVec3d(origin);
  vp.CopyTo(p);
}

void delfem2::opengl::CRender2Tex_DrawOldGL::BoundingBox3
 (double* pmin,
  double* pmax) const
{
  for(unsigned int ix=0;ix<nResX;++ix){
    for(unsigned int iy=0;iy<nResY;++iy){
      CVec3d vp;
      {
        const CVec3d& dx = x_axis;
        const CVec3d& dz = z_axis;
        const CVec3d& dy = Cross(dz,dx);
        double lz = aZ[iy*nResX+ix];
        double lx = (ix+0.5)*lengrid;
        double ly = (iy+0.5)*lengrid;
        vp = lx*dx+ly*dy+lz*dz + CVec3d(origin);
        if( -lz > z_range*0.99 ) continue;
      }
      // ----------
      if( pmin[0] > pmax[0] ){
        pmin[0] = pmax[0] = vp.x();
        pmin[1] = pmax[1] = vp.y();
        pmin[2] = pmax[2] = vp.z();
        continue;
      }
      const double x0 = vp.x();
      if( x0 < pmin[0] ){ pmin[0] = x0; }
      if( x0 > pmax[0] ){ pmax[0] = x0; }
      const double y0 = vp.y();
      if( y0 < pmin[1] ){ pmin[1] = y0; }
      if( y0 > pmax[1] ){ pmax[1] = y0; }
      const double z0 = vp.z();
      if( z0 < pmin[2] ){ pmin[2] = z0; }
      if( z0 > pmax[2] ){ pmax[2] = z0; }
    }
  }
}


// --------------------------------------------------------

void delfem2::opengl::CRender2Tex_DrawOldGL_BOX::Initialize
 (unsigned int nresX_,
  unsigned int nresY_,
  unsigned int nresZ_,
  double elen_)
{
  aSampler.resize(6);
  aSampler[0].SetTextureProperty(nresY_, nresZ_, true); // +x
  aSampler[0].SetCoord(elen_, elen_*nresX_,
                       CVec3d(+0.5*elen_*nresX_,-0.5*elen_*nresY_,-0.5*elen_*nresZ_).stlvec(),
                       CVec3d(+1,  0, 0).stlvec(),
                       CVec3d( 0, +1, 0).stlvec() );
  aSampler[0].SetPointColor(1.0, 0.0, 0.0);
  //
  aSampler[1].SetTextureProperty(nresY_, nresZ_, true); // -x
  aSampler[1].SetCoord(elen_, elen_*nresX_,
                       CVec3d(-0.5*elen_*nresX_,-0.5*elen_*nresY_,+0.5*elen_*nresZ_).stlvec(),
                       CVec3d(-1,  0, 0).stlvec(),
                       CVec3d( 0, +1, 0).stlvec() );
  aSampler[1].SetPointColor(1.0, 0.5, 0.5);
  //
  aSampler[2].SetTextureProperty(nresX_, nresZ_, true); // +y
  aSampler[2].SetCoord(elen_, elen_*nresY_,
                       CVec3d(-0.5*elen_*nresX_,+0.5*elen_*nresY_,+0.5*elen_*nresZ_).stlvec(),
                       CVec3d(0,+1,0).stlvec(),
                       CVec3d(1,+0,0).stlvec() );
  aSampler[2].SetPointColor(0.0, 1.0, 0.0);
  //
  aSampler[3].SetTextureProperty(nresX_, nresZ_, true); // -y
  aSampler[3].SetCoord(elen_, elen_*nresY_,
                       CVec3d(-0.5*elen_*nresX_,-0.5*elen_*nresY_,-0.5*elen_*nresZ_).stlvec(),
                       CVec3d(0,-1,0).stlvec(),
                       CVec3d(1,+0,0).stlvec() );
  aSampler[3].SetPointColor(0.5, 1.0, 0.5);
  //
  aSampler[4].SetTextureProperty(nresX_, nresY_, true);
  aSampler[4].SetCoord(elen_, elen_*nresZ_,
                       CVec3d(-0.5*elen_*nresX_,-0.5*elen_*nresY_,+0.5*elen_*nresZ_).stlvec(),
                       CVec3d(0,0,+1).stlvec(),
                       CVec3d(1,0,0).stlvec() );
  aSampler[4].SetPointColor(0.0, 0.0, 1.0);
  //
  aSampler[5].SetTextureProperty(nresX_, nresY_, true);
  aSampler[5].SetCoord(elen_, elen_*nresZ_,
                       CVec3d(-0.5*elen_*nresX_,+0.5*elen_*nresY_,-0.5*elen_*nresZ_).stlvec(),
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


void delfem2::opengl::CarveVoxelByDepth
 (std::vector<int>& aVal,
  const CRender2Tex_DrawOldGL_BOX& sampler_box)
{
  const unsigned int nx = sampler_box.nDivX();
  const unsigned int ny = sampler_box.nDivY();
  const unsigned int nz = sampler_box.nDivZ();
  const double el = sampler_box.edgeLen();
  // ------
  aVal.assign(nz*ny*nx,1);
  for(unsigned int iy=0;iy<ny;++iy){
    for(unsigned int iz=0;iz<nz;++iz){
      double d0 = sampler_box.aSampler[0].aZ[iz*ny+iy];
      const unsigned int nd = floor(-d0/el+1.0e-5);
      for(unsigned int id=0;id<nd;id++){
        const unsigned int ix0 = nx-1-id;
        const unsigned int iy0 = iy;
        const unsigned int iz0 = iz;
        aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
      }
    }
  }
  for(unsigned int iy=0;iy<ny;++iy){
    for(unsigned int iz=0;iz<nz;++iz){
      double d0 = sampler_box.aSampler[1].aZ[iz*ny+iy];
      const unsigned int nd = floor(-d0/el+1.0e-5);
      for(unsigned int id=0;id<nd;id++){
        const unsigned int ix0 = id;
        const unsigned int iy0 = iy;
        const unsigned int iz0 = nz-1-iz;
        aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
      }
    }
  }
  for(unsigned int ix=0;ix<nx;++ix){
    for(unsigned int iz=0;iz<nz;++iz){
      double d0 = sampler_box.aSampler[2].aZ[iz*nx+ix];
      const unsigned int nd = floor(-d0/el+1.0e-5);
      for(unsigned int id=0;id<nd;id++){
        const unsigned int ix0 = ix;
        const unsigned int iy0 = ny-1-id;
        const unsigned int iz0 = nz-1-iz;
        aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
      }
    }
  }
  for(unsigned int ix=0;ix<nx;++ix){
    for(unsigned int iz=0;iz<nz;++iz){
      double d0 = sampler_box.aSampler[3].aZ[iz*nx+ix];
      const unsigned int nd = floor(-d0/el+1.0e-5);
      for(unsigned int id=0;id<nd;id++){
        const unsigned int ix0 = ix;
        const unsigned int iy0 = id;
        const unsigned int iz0 = iz;
        aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
      }
    }
  }
  for(unsigned int ix=0;ix<nx;++ix){
    for(unsigned int iy=0;iy<ny;++iy){
      double d0 = sampler_box.aSampler[4].aZ[iy*nx+ix];
      const unsigned int nd = floor(-d0/el+1.0e-5);
      for(unsigned int id=0;id<nd;id++){
        const unsigned int ix0 = ix;
        const unsigned int iy0 = iy;
        const unsigned int iz0 = nz-1-id;
        aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
      }
    }
  }
  for(unsigned int ix=0;ix<nx;++ix){
    for(unsigned int iy=0;iy<ny;++iy){
      double d0 = sampler_box.aSampler[5].aZ[iy*nx+ix];
      const unsigned int nd = floor(-d0/el+1.0e-5);
      for(unsigned int id=0;id<nd;id++){
        const unsigned int ix0 = ix;
        const unsigned int iy0 = ny-1-iy;
        const unsigned int iz0 = id;
        aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
      }
    }
  }
}
