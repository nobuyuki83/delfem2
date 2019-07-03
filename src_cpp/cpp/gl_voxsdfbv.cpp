/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#if defined(__APPLE__) && defined(__MACH__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else // linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif


#include "delfem2/vec3.h"

#include "delfem2/gl_voxsdfbv.h"


static void myGlNormal(const CVector3& n){ ::glNormal3d(n.x,n.y,n.z); }
static void myGlVertex(const CVector3& v){ ::glVertex3d(v.x,v.y,v.z); }
static void myGlColorDiffuse(float r, float g, float b, float a){
  ::glColor4d(r, g, b, a );
  float c[4] = {r, g, b, a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

/////////////////////////////

void Draw(const CSignedDistanceField3D_Plane& sdf)
{
  const double* origin_ = sdf.origin_;
  const double* normal_ = sdf.normal_;
  //////
  const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(1,0,0);
  double d1[3],d2[3];
  GetVertical2Vector3D(normal_,d1,d2);
  ::glBegin(GL_LINES);
  for(int i=-10;i<11;i++){
    for(int j=-10;j<11;j++){
      ::glVertex3d(origin_[0]+d1[0]*0.08*(i+0)+d2[0]*0.08*(j+0),
                   origin_[1]+d1[1]*0.08*(i+0)+d2[1]*0.08*(j+0),
                   origin_[2]+d1[2]*0.08*(i+0)+d2[2]*0.08*(j+0));
      ::glVertex3d(origin_[0]+d1[0]*0.08*(i+1)+d2[0]*0.08*(j+0),
                   origin_[1]+d1[1]*0.08*(i+1)+d2[1]*0.08*(j+0),
                   origin_[2]+d1[2]*0.08*(i+1)+d2[2]*0.08*(j+0));
      
      ::glVertex3d(origin_[0]+d1[0]*0.08*(i+0)+d2[0]*0.08*(j+0),
                   origin_[1]+d1[1]*0.08*(i+0)+d2[1]*0.08*(j+0),
                   origin_[2]+d1[2]*0.08*(i+0)+d2[2]*0.08*(j+0));
      ::glVertex3d(origin_[0]+d1[0]*0.08*(i+0)+d2[0]*0.08*(j+1),
                   origin_[1]+d1[1]*0.08*(i+0)+d2[1]*0.08*(j+1),
                   origin_[2]+d1[2]*0.08*(i+0)+d2[2]*0.08*(j+1));
    }
  }
  ::glEnd();
  /*  const double para[3] = { normal_[1], -normal_[0], 0 };
   ::glColor3d(0,0,0);
   ::glBegin(GL_QUADS);
   ::glVertex3d(origin_[0]+para[0]-normal_[0], origin_[1]+para[1]-normal_[1], origin_[2]+para[2]-normal_[2] );
   ::glVertex3d(origin_[0]-para[0]-normal_[0], origin_[1]-para[1]-normal_[1], origin_[2]-para[2]-normal_[2] );
   ::glVertex3d(origin_[0]-para[0], origin_[1]-para[1], origin_[2]-para[2] );
   ::glVertex3d(origin_[0]+para[0], origin_[1]+para[1], origin_[2]+para[2] );
   ::glEnd();*/
  if(is_lighting){ glEnable(GL_LIGHTING); }
}



void DrawSDF_Sphere(const CSignedDistanceField3D_Sphere& sdf)
{
  const double* cent_ = sdf.cent_;
  const double radius_ = sdf.radius_;
  ////////
  const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  const bool is_texture = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_TEXTURE_2D);
  ////
  ::glLineWidth(1);
  ::glColor3d(1,0,0);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glTranslated(cent_[0],cent_[1],cent_[2]);
  
  const unsigned int nlg = 32;
  const unsigned int nlt = 18;
  const double rlg = 6.28/nlg;  // longtitude
  //  const double rlt = 6.28/nlt;  // latitude
  //  const double divlg = 6.28/ndiv_lg;
  //  const double divlt = 6.28/ndiv_lt;
  const unsigned int ndiv = 32;
  const double rdiv = 6.28/ndiv;
  for(unsigned int ilg=0;ilg<nlg;ilg++){
    ::glBegin(GL_LINE_LOOP);
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      ::glVertex3d(radius_*cos(idiv*rdiv)*cos(ilg*rlg),
                   radius_*cos(idiv*rdiv)*sin(ilg*rlg),
                   radius_*sin(idiv*rdiv) );
    }
    ::glEnd();
  }
  for(unsigned int ilt=0;ilt<nlt;ilt++){
    const double d = ((double)ilt/nlt-0.5)*radius_*2.0;
    const double r0 = sqrt(radius_*radius_-d*d);
    ::glBegin(GL_LINE_LOOP);
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      ::glVertex3d(r0*cos(idiv*rdiv),
                   r0*sin(idiv*rdiv),
                   d);
    }
    ::glEnd();
  }
  //   ::glutWireSphere(radius_,32,32);
  ::glPopMatrix();
  ////
  if(is_lighting){ glEnable(GL_LIGHTING); }
  if(is_texture ){ glEnable(GL_TEXTURE_2D); }
}


void DrawFace(const CSignedDistanceField3D_Cylinder& sdf)
{
  const double* dir_ = sdf.dir_;
  const double radius_ = sdf.radius_;
  const double* cent_ = sdf.cent_;
  
  ////////////
  //  const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  const bool is_texture = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_TEXTURE_2D);
  ////
  ::glLineWidth(1);
  ::glColor3d(1,0,0);
  float brown[3] = {0.7,0.20,0.15};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, brown);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  
  double h[3],v[3];
  {
    double x[3] = { 1,0,0 };
    Cross3D(h, dir_, x);
    if( SquareLength3D(h) > 1.0e-5 ){
      double invlenh=1.0/Length3D(h);
      h[0]*=invlenh;
      h[1]*=invlenh;
      h[2]*=invlenh;
      Cross3D(v, h, dir_);
    }
    else{
      double y[3] = { 0,1,0 };
      Cross3D(h, dir_, y);
      double invlenh=1.0/Length3D(h);
      h[0]*=invlenh;
      h[1]*=invlenh;
      h[2]*=invlenh;
      Cross3D(v, h, dir_);
    }
  }
  
  const unsigned int nl = 20;
  const unsigned int nr = 32;
  const double rr = 6.28/nr;
  
  for(unsigned int ir=0;ir<nr;ir++){
    double rs0 = +radius_*sin(ir*rr);
    double rc0 = +radius_*cos(ir*rr);
    double rs1 = +radius_*sin(ir*rr+rr);
    double rc1 = +radius_*cos(ir*rr+rr);
    
    double v0[3] = {
      cent_[0]+rs0*h[0]+rc0*v[0],
      cent_[1]+rs0*h[1]+rc0*v[1],
      cent_[2]+rs0*h[2]+rc0*v[2] };
    double v1[3] = {
      cent_[0]+rs1*h[0]+rc1*v[0],
      cent_[1]+rs1*h[1]+rc1*v[1],
      cent_[2]+rs1*h[2]+rc1*v[2] };
    ::glBegin(GL_TRIANGLES);
    ::glNormal3d(+dir_[0],+dir_[1],+dir_[2]);
    ::glVertex3d(cent_[0]+dir_[0], cent_[1]+dir_[1], cent_[2]+dir_[2]);
    ::glVertex3d(dir_[0]+v0[0], dir_[1]+v0[1], dir_[2]+v0[2]);
    ::glVertex3d(dir_[0]+v1[0], dir_[1]+v1[1], dir_[2]+v1[2]);
    /////
    ::glNormal3d(-dir_[0],-dir_[1],-dir_[2]);
    ::glVertex3d(cent_[0]-dir_[0], cent_[1]-dir_[1], cent_[2]-dir_[2]);
    ::glVertex3d(-dir_[0]+v0[0],-dir_[1]+v0[1],-dir_[2]+v0[2]);
    ::glVertex3d(-dir_[0]+v1[0],-dir_[1]+v1[1],-dir_[2]+v1[2]);
    ::glEnd();
    /////
    double s0 = sin(ir*rr);
    double c0 = cos(ir*rr);
    double n0[3] = {
      s0*h[0]+c0*v[0],
      s0*h[1]+c0*v[1],
      s0*h[2]+c0*v[2] };
    /////
    double s1 = sin(ir*rr+rr);
    double c1 = cos(ir*rr+rr);
    double n1[3] = {
      s1*h[0]+c1*v[0],
      s1*h[1]+c1*v[1],
      s1*h[2]+c1*v[2] };
    
    ::glBegin(GL_QUADS);
    for(unsigned int il=0;il<nl;il++){
      double l0 = -1+(2.0/nl)*il;
      double l1 = -1+(2.0/nl)*(il+1);
      ::glNormal3d(+n0[0],+n0[1],+n0[2]);
      ::glVertex3d(l0*dir_[0]+v0[0], l0*dir_[1]+v0[1], l0*dir_[2]+v0[2]);
      ::glVertex3d(l1*dir_[0]+v0[0], l1*dir_[1]+v0[1], l1*dir_[2]+v0[2]);
      ::glNormal3d(+n1[0],+n1[1],+n1[2]);
      ::glVertex3d(l1*dir_[0]+v1[0], l1*dir_[1]+v1[1], l1*dir_[2]+v1[2]);
      ::glVertex3d(l0*dir_[0]+v1[0], l0*dir_[1]+v1[1], l0*dir_[2]+v1[2]);
    }
    ::glEnd();
  }
  ::glPopMatrix();
  ////
  if(is_texture ){ glEnable(GL_TEXTURE_2D); }
}


void Draw(const CSignedDistanceField3D_Cylinder& sdf)
{
  const double* dir_ = sdf.dir_;
  const double radius_ = sdf.radius_;
  const double* cent_ = sdf.cent_;
  
  //////
  const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  const bool is_texture = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_TEXTURE_2D);
  ////
  ::glLineWidth(1);
  ::glColor3d(1,0,0);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  
  double h[3],v[3];
  {
    double x[3] = { 1,0,0 };
    Cross3D(h, dir_, x);
    if( SquareLength3D(h) > 1.0e-5 ){
      double invlenh=1.0/Length3D(h);
      h[0]*=invlenh;
      h[1]*=invlenh;
      h[2]*=invlenh;
      Cross3D(v, h, dir_);
    }
    else{
      double y[3] = { 0,1,0 };
      Cross3D(h, dir_, y);
      double invlenh=1.0/Length3D(h);
      h[0]*=invlenh;
      h[1]*=invlenh;
      h[2]*=invlenh;
      Cross3D(v, h, dir_);
    }
  }
  
  const unsigned int nr = 32;
  const double rr = 6.28/nr;
  ::glBegin(GL_LINES);
  for(unsigned int ir=0;ir<nr;ir++){
    double rs = +radius_*sin(ir*rr);
    double rc = +radius_*cos(ir*rr);
    ::glVertex3d(cent_[0]+dir_[0],
                 cent_[1]+dir_[1],
                 cent_[2]+dir_[2]);
    ::glVertex3d(cent_[0]+dir_[0]+rs*h[0]+rc*v[0],
                 cent_[1]+dir_[1]+rs*h[1]+rc*v[1],
                 cent_[2]+dir_[2]+rs*h[2]+rc*v[2]);
    ::glVertex3d(cent_[0]-dir_[0],
                 cent_[1]-dir_[1],
                 cent_[2]-dir_[2]);
    ::glVertex3d(cent_[0]-dir_[0]+rs*h[0]+rc*v[0],
                 cent_[1]-dir_[1]+rs*h[1]+rc*v[1],
                 cent_[2]-dir_[2]+rs*h[2]+rc*v[2]);
    ::glVertex3d(cent_[0]+dir_[0]+rs*h[0]+rc*v[0],
                 cent_[1]+dir_[1]+rs*h[1]+rc*v[1],
                 cent_[2]+dir_[2]+rs*h[2]+rc*v[2]);
    ::glVertex3d(cent_[0]-dir_[0]+rs*h[0]+rc*v[0],
                 cent_[1]-dir_[1]+rs*h[1]+rc*v[1],
                 cent_[2]-dir_[2]+rs*h[2]+rc*v[2]);
  }
  ::glEnd();
  ::glPopMatrix();
  ////
  if(is_lighting){ glEnable(GL_LIGHTING); }
  if(is_texture ){ glEnable(GL_TEXTURE_2D); }
}


void Draw(const CSignedDistanceField3D_Torus& sdf)
{
  const double radius_ = sdf.radius_;
  const double radius_tube_ = sdf.radius_tube_;
  const double* cent_ = sdf.cent_;
  ////////
  const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  const bool is_texture = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_TEXTURE_2D);
  ////
  ::glColor3d(1,0,0);
  ::glPushMatrix();
  ::glTranslated(cent_[0],cent_[1],cent_[2]);
  
  const unsigned int nlg = 32;
  const unsigned int nlt = 18;
  const double rlg = 6.28/nlg;  // longtitude
  const double rlt = 6.28/nlt;  // latitude
  const unsigned int ndiv = 32;
  const double rdiv = 6.28/ndiv;
  for(unsigned int ilg=0;ilg<nlg;ilg++){
    ::glBegin(GL_LINE_LOOP);
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      ::glVertex3d(( radius_ + radius_tube_*cos(idiv*rdiv) )*sin(ilg*rlg),
                   ( radius_ + radius_tube_*cos(idiv*rdiv) )*cos(ilg*rlg),
                   radius_tube_*sin(idiv*rdiv) );
    }
    ::glEnd();
  }
  for(unsigned int ilt=0;ilt<nlt;ilt++){
    double d  = radius_tube_*sin(ilt*rlt);
    double r0 = radius_tube_*cos(ilt*rlt) + radius_;
    ::glBegin(GL_LINE_LOOP);
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      ::glVertex3d(r0*cos(idiv*rdiv),
                   r0*sin(idiv*rdiv),
                   d);
    }
    ::glEnd();
  }
  ::glPopMatrix();
  ////
  if(is_lighting){ ::glEnable(GL_LIGHTING); }
  if(is_texture) { ::glEnable(GL_TEXTURE_2D); }
}


void Draw(const CSignedDistanceField3D_Mesh& sdf)
{
  const unsigned int* aTri_ = sdf.aTri_;
  const double* pXYZs_ = sdf.pXYZs_;
  const int ntri_ = sdf.ntri_;
  ////////
  ::glColor3d(1,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itri=0;itri<ntri_;itri++){
    int i1 = aTri_[itri*3+0];
    int i2 = aTri_[itri*3+1];
    int i3 = aTri_[itri*3+2];
    ::glVertex3dv(pXYZs_+i1*3);
    ::glVertex3dv(pXYZs_+i2*3);
    ::glVertex3dv(pXYZs_+i2*3);
    ::glVertex3dv(pXYZs_+i3*3);
    ::glVertex3dv(pXYZs_+i3*3);
    ::glVertex3dv(pXYZs_+i1*3);
  }
  ::glEnd();
}


////////////////////////////////////////////////////////



const unsigned int noelElemFace_Vox[6][4] = {
  { 0, 4, 6, 2 }, // -x
  { 1, 3, 7, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 2, 6, 7, 3 }, // +y
  { 0, 2, 3, 1 }, // -z
  { 4, 5, 7, 6 } }; // +z

const CVector3 normalHexFace[6] = {
  CVector3(-1, 0, 0),
  CVector3(+1, 0, 0),
  CVector3( 0,-1, 0),
  CVector3( 0,+1, 0),
  CVector3( 0, 0,-1),
  CVector3( 0, 0,+1)
};

void Draw_CubeGrid
(bool is_picked, int iface_picked,
 double elen, const CVector3& org,
 const CCubeGrid& cube)
{
  if( !cube.is_active ) return;
  int ih = cube.ivx;
  int jh = cube.ivy;
  int kh = cube.ivz;
  CVector3 aP[8] = {
    org + elen*CVector3(ih+0,jh+0,kh+0),
    org + elen*CVector3(ih+1,jh+0,kh+0),
    org + elen*CVector3(ih+0,jh+1,kh+0),
    org + elen*CVector3(ih+1,jh+1,kh+0),
    org + elen*CVector3(ih+0,jh+0,kh+1),
    org + elen*CVector3(ih+1,jh+0,kh+1),
    org + elen*CVector3(ih+0,jh+1,kh+1),
    org + elen*CVector3(ih+1,jh+1,kh+1) };
  ::glEnable(GL_LIGHTING);
  ::glBegin(GL_QUADS);
  for(int iface=0;iface<6;++iface){
    if( is_picked && iface_picked == iface ){
      ::myGlColorDiffuse(1,1,0,1);
    }
    else{
      ::myGlColorDiffuse(0.8,0.8,0.8,1);
    }
    myGlNormal(normalHexFace[iface]);
    myGlVertex(aP[noelElemFace_Vox[iface][0]]);
    myGlVertex(aP[noelElemFace_Vox[iface][1]]);
    myGlVertex(aP[noelElemFace_Vox[iface][2]]);
    myGlVertex(aP[noelElemFace_Vox[iface][3]]);
  }
  ::glEnd();
  //////////////////////////////
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



void Draw(const CBV3D_AABB& aabb)
{
  double x_min = aabb.x_min;
  double x_max = aabb.x_max;
  double y_min = aabb.y_min;
  double y_max = aabb.y_max;
  double z_min = aabb.z_min;
  double z_max = aabb.z_max;
  /////
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


