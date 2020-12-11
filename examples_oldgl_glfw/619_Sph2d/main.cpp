/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @file implementation based on "Müller et al., Particle-based fluid simulation for interactive applications. SCA 2003"
 */

// ----------------------
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "delfem2/vec3.h"

// -----------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

#ifndef M_PI
#  define M_PI 3.1415926535
#endif

namespace dfm2 = delfem2;

// -----------------

struct SParticle
{
  dfm2::CVec3d r;  // position
  dfm2::CVec3d v;  // velocity
  double rho;   // mass-density
  double p;
  dfm2::CVec3d f;
};

void SPH_DensityPressure
 (std::vector<SParticle>& ps,
  double H,
  double SPH_SIMSCALE,
  double SPH_PMASS,
  double SPH_RESTDENSITY,
  double SPH_INTSTIFF)
{
  const double Poly6Kern = 315.0 / ( 64.0 * M_PI * pow( H, 9 ) );
//  const double SpikyKern = -45.0 / ( M_PI * pow( H, 6 ) );
//  const double LapKern   = 45.0 / ( M_PI * pow( H, 6 ) );
  for(unsigned int ips=0;ips<ps.size();ips++){
    SParticle& ps0 = ps[ips];
    double sum = 0.0;
    for(unsigned int jps=0;jps<ps.size();jps++){
      if ( ips == jps ) continue;
      const SParticle& ps1 = ps[jps];
      double dr[3] = {
        (ps0.r[0]-ps1.r[0])*SPH_SIMSCALE,
        (ps0.r[1]-ps1.r[1])*SPH_SIMSCALE,
        (ps0.r[2]-ps1.r[2])*SPH_SIMSCALE };
      const double sqr = dfm2::SquareLength3( dr );
      if ( H*H < sqr ) continue;
      const double c = H*H - sqr;
      sum += c * c * c;
    }
    ps0.rho = sum * SPH_PMASS * Poly6Kern;
    //    std::cout << "rho " << ps0.rho << std::endl;
    ps0.p   = ( ps0.rho - SPH_RESTDENSITY ) * SPH_INTSTIFF;
    //    std::cout << "p   " << ps0.p << std::endl;
    ps0.rho = 1.0 / ps0.rho;  // take inverse for later calculation
  }
}

void SPH_Force
(std::vector<SParticle>& ps,
 double H,
 double SPH_SIMSCALE,
 double SPH_VISC)
{
//  const double Poly6Kern = 315.0 / ( 64.0 * M_PI * pow( H, 9 ) );
  const double SpikyKern = -45.0 / ( M_PI * pow( H, 6 ) );
  const double LapKern   = 45.0 / ( M_PI * pow( H, 6 ) );
  for(unsigned int ips=0;ips<ps.size();ips++){
    SParticle& ps0 = ps[ips];
    double force[3] = { 0,0,0 };
    for(unsigned int jps=0;jps<ps.size();jps++){
      if ( ips == jps ) continue;
      const SParticle& ps1 = ps[jps];
      double dr[3] = {
        ( ps0.r[0] - ps1.r[0] ) * SPH_SIMSCALE,
        ( ps0.r[1] - ps1.r[1] ) * SPH_SIMSCALE,
        ( ps0.r[2] - ps1.r[2] ) * SPH_SIMSCALE };
      const double r = dfm2::Length3( dr );
      if ( H < r ) continue;
      const double c = H - r;
      const double pterm = -0.5 * c * SpikyKern * ( ps0.p + ps1.p ) / r;
      const double vterm = LapKern * SPH_VISC;
      double fcurr[3] = {
        pterm * dr[0] + vterm * ( ps1.v[0] - ps0.v[0] ),
        pterm * dr[1] + vterm * ( ps1.v[1] - ps0.v[1] ),
        pterm * dr[2] + vterm * ( ps1.v[2] - ps0.v[2] ) };
      fcurr[0] *= c * ps0.rho * ps1.rho;
      fcurr[1] *= c * ps0.rho * ps1.rho;
      fcurr[2] *= c * ps0.rho * ps1.rho;
      force[0] += fcurr[0];
      force[1] += fcurr[1];
      force[2] += fcurr[2];
    }
    //    std::cout << "force " << force[0] << " " << force[1] << " " << force[2] << std::endl;
    ps0.f[0] = force[0];
    ps0.f[1] = force[1];
    ps0.f[2] = force[2];
  }
}

void SPH_UpdatePosition
 (std::vector<SParticle>& ps,
  double SPH_PMASS,
  double SPH_LIMIT,
  double SPH_SIMSCALE,
  double SPH_EXTSTIFF,
  double SPH_EXTDAMP,
  double DT,
  double EPSILON,
  double MIN[3],
  double MAX[3],
  double RADIUS)
{
  const double g[3] = { 0.0, -9.8, 0.0};
  for(unsigned int ips=0;ips<ps.size();ips++){
    SParticle& p = ps[ips];
    double accel[3] = {
      p.f[0] * SPH_PMASS,
      p.f[1] * SPH_PMASS,
      p.f[2] * SPH_PMASS };
    const double sq_speed = dfm2::SquareLength3( accel );
    if ( sq_speed > SPH_LIMIT*SPH_LIMIT ) {
      accel[0] *= SPH_LIMIT / sqrt(sq_speed);
      accel[1] *= SPH_LIMIT / sqrt(sq_speed);
      accel[2] *= SPH_LIMIT / sqrt(sq_speed);
    }
    for(unsigned int idim=0;idim<3;idim++){
      {
        const double diff = 2.0 * RADIUS - ( p.r[idim] - MIN[idim] ) * SPH_SIMSCALE;
        if ( diff > EPSILON )
        {
          double norm[3] = { 0.0, 0.0, 0.0 }; norm[idim] = +1.0;
          const double adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * dfm2::Dot3( norm, p.v.p );
          accel[0] += adj * norm[0];
          accel[1] += adj * norm[1];
          accel[2] += adj * norm[2];
        }
      }
      {
        const double diff = 2.0 * RADIUS - ( MAX[idim] - p.r[idim] ) * SPH_SIMSCALE;
        if ( diff > EPSILON )
        {
          double norm[3] = { 0.0, 0.0, 0.0 }; norm[idim] = -1.0;
          const double adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * dfm2::Dot3( norm, p.v.p );
          accel[0] += adj * norm[0];
          accel[1] += adj * norm[1];
          accel[2] += adj * norm[2];
        }
      }
    }
    
    accel[0] += g[0];
    accel[1] += g[1];
    accel[2] += g[2];
    
    p.v[0] += accel[0] * DT;
    p.v[1] += accel[1] * DT;
    p.v[2] += accel[2] * DT;
    
    p.r[0] += p.v[0] * DT / SPH_SIMSCALE;
    p.r[1] += p.v[1] * DT / SPH_SIMSCALE;
    p.r[2] += p.v[2] * DT / SPH_SIMSCALE;
    
    // stick on x-y plane
    p.r[2] = 0.0;
  }
}

// --------------------------------

double H;
double DT;
double RADIUS;
double EPSILON;
double MIN[3];
double MAX[3];
double SPH_PMASS;
double SPH_INTSTIFF;
double SPH_EXTSTIFF;
double SPH_EXTDAMP;
double SPH_PDIST;
double SPH_SIMSCALE;
double SPH_RESTDENSITY;
double SPH_VISC;
double SPH_LIMIT;
std::vector<SParticle> ps;

void init( std::vector<SParticle>& ps )
{
  H         = 0.01; // m
  DT        = 0.001;
  RADIUS    = 0.001; // m
  EPSILON   = 0.00001;
  // size of water column
  // size of box
  MIN[0] =  0.0; MIN[1] =  0.0; MIN[2] = -10.0;
  MAX[0] = 20.0; MAX[1] = 20.0; MAX[2] =  10.0;
  SPH_PMASS = 0.00020543; // kg
  SPH_INTSTIFF = 1.00;  //
  SPH_EXTSTIFF = 50000.0; // penalty coefficient of the wall repulsion force
  SPH_EXTDAMP  = 256.0; // damping of the wall repulsion force
  SPH_RESTDENSITY = 600.0; // kg / m^3
  SPH_PDIST = pow( SPH_PMASS / SPH_RESTDENSITY, 1/3.0 );
  SPH_SIMSCALE = 0.008; // unit size
  SPH_VISC = 0.2; // pascal-second (Pa.s) = 1 kg m^-1 s^-1
  SPH_LIMIT = 200.0;
  
  double d;
  d = SPH_PDIST * 0.87 / SPH_SIMSCALE;
  
  ps.clear();
  {
    const double INITMIN[3] = { 0.0, 0.0, 0.0};
    const double INITMAX[3] = { 10.0, 20.0, 0.0};
    for ( double z = INITMIN[2]; z <= INITMAX[2]; z += d ) {
    for ( double y = INITMIN[1]; y <= INITMAX[1]; y += d ) {
    for ( double x = INITMIN[0]; x <= INITMAX[0]; x += d ) {
      SParticle p;
      p.r[0] = x; p.r[1] = y; p.r[2] = z;
      p.r[0] += -0.05 + rand() * 0.1 / RAND_MAX;
      p.r[1] += -0.05 + rand() * 0.1 / RAND_MAX;
      p.r[2] += -0.05 + rand() * 0.1 / RAND_MAX;
      p.v[0] = 0.0; p.v[1] = 0.0; p.v[2] = 0.0;
      ps.push_back( p );
    }
    }
    }
  }
}
      
void myGlutDisplay(void)
{  
  ::glClearColor(1,1,1,1);
  
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);
  
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
  ::glDisable(GL_LIGHTING);
  
  ::glColor3d(0,0,0);
  ::glPointSize(5);
  ::glBegin(GL_POINTS);
  for(unsigned int ips=0;ips<ps.size();ips++){
    ::glVertex2dv(ps[ips].r.p);
  }
  ::glEnd();
  ::glBegin(GL_LINES);
  ::glVertex3d( 0, 0,0);
  ::glVertex3d( 0,20,0);
  ::glVertex3d( 0,20,0);
  ::glVertex3d(20,20,0);  
  ::glVertex3d(20,20,0);  
  ::glVertex3d(20, 0,0);
  ::glVertex3d(20, 0,0);  
  ::glVertex3d( 0, 0,0);    
  ::glEnd();
}
  
int main(int argc, char *argv[])
{
  init( ps );
  std::cout << "particle size : " << ps.size() << std::endl;  
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 25;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  
  while (true){
    {
      SPH_DensityPressure( ps,
                          H,SPH_SIMSCALE,SPH_PMASS,SPH_RESTDENSITY,SPH_INTSTIFF);
      SPH_Force(ps,
                H,SPH_SIMSCALE,SPH_VISC);
      SPH_UpdatePosition( ps,
                         SPH_PMASS,SPH_LIMIT,SPH_SIMSCALE,SPH_EXTSTIFF,SPH_EXTDAMP,
                         DT,EPSILON,MIN,MAX,RADIUS);
    }
    // -----
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto CLOSE; }
  }
  
CLOSE:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
