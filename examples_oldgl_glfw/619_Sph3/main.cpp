/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @file implementation based on "Müller et al., Particle-based fluid simulation for interactive applications. SCA 2003"
 */

#include <cmath>
#include <iostream>
#include <vector>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/vec3.h"
#include "delfem2/srchgrid.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"


#ifndef M_PI
  #define M_PI 3.1415926535
#endif

namespace dfm2 = delfem2;

// -----------------

struct SParticle
{
  dfm2::CVec3d r;  //! position
  dfm2::CVec3d v;  //! velocity
  double rho;   //! mass-density
  double p; //! pressure
  dfm2::CVec3d f; //! force
};

double SPH_CubicSquareDistance(
    const SParticle& ps0,
    const SParticle& ps1,
    double sph_radcutoff)
{
  const double dr[3] = {
      ps0.r[0]-ps1.r[0],
      ps0.r[1]-ps1.r[1],
      ps0.r[2]-ps1.r[2] };
  const double sqr = dfm2::SquareLength3( dr );
  if (sph_radcutoff * sph_radcutoff < sqr ) { return 0.0; }
  const double c = sph_radcutoff * sph_radcutoff - sqr;
  return c * c * c;
}

/**
 * update density&pressure of the particle
 * @param ps
 * @param H
 * @param SPH_SIMSCALE
 * @param SPH_PMASS
 * @param SPH_RESTDENSITY
 * @param SPH_INTSTIFF
 */
void SPH_DensityPressure(
    std::vector<SParticle>& ps,
    double H,
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
      sum += SPH_CubicSquareDistance(ps0, ps1, H);
    }
    ps0.rho = sum * SPH_PMASS * Poly6Kern;
    ps0.p   = ( ps0.rho - SPH_RESTDENSITY ) * SPH_INTSTIFF;
    ps0.rho = 1.0 / ps0.rho;  // take inverse for later calculation
  }
}

template <class Search>
void SPH_DensityPressureHash(
    std::vector<SParticle>& ps,
    double sph_radcutoff,
    double psh_pmass,
    double sph_restdensity,
    double sph_intstiff,
    const Search& sg)
{
  const double Poly6Kern = 315.0 / ( 64.0 * M_PI * pow(sph_radcutoff, 9 ) );
//  const double SpikyKern = -45.0 / ( M_PI * pow( H, 6 ) );
//  const double LapKern   = 45.0 / ( M_PI * pow( H, 6 ) );
  std::vector<unsigned int> aIP;
  for(unsigned int ips=0;ips<ps.size();ips++){
    SParticle& ps0 = ps[ips];
    double sum = 0.0;
    sg.GetOneRingNeighbor(aIP,ps0.r.p);
    for(unsigned int jps : aIP) {
      if( ips == jps ){ continue; }
      sum += SPH_CubicSquareDistance(ps0, ps[jps], sph_radcutoff);
    }
    ps0.rho = sum * psh_pmass * Poly6Kern;
    ps0.p   = ( ps0.rho - sph_restdensity ) * sph_intstiff;
    ps0.rho = 1.0 / ps0.rho;  // take inverse for later calculation
  }
}

void Force(
    double force[3],
    const SParticle& ps0,
    const SParticle& ps1,
    double sph_radcutoff,
    double SpikyKern,
    double vterm)
{
  double dr[3] = {
      ps0.r[0] - ps1.r[0],
      ps0.r[1] - ps1.r[1],
      ps0.r[2] - ps1.r[2] };
  const double r = dfm2::Length3( dr );
  if ( sph_radcutoff < r ) return;
  const double c = sph_radcutoff - r;
  const double pterm = -0.5 * c * SpikyKern * ( ps0.p + ps1.p ) / r;
//  const double vterm = LapKern * SPH_VISC;
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

void SPH_AccumulateForce(
    std::vector<SParticle>& ps,
    double sph_radcutoff,
    double SPH_VISC)
{
//  const double Poly6Kern = 315.0 / ( 64.0 * M_PI * pow( H, 9 ) );
  const double SpikyKern = -45.0 / ( M_PI * pow( sph_radcutoff, 6 ) );
  const double LapKern   = 45.0 / ( M_PI * pow( sph_radcutoff, 6 ) );
  for(unsigned int ips=0;ips<ps.size();ips++){
    SParticle& ps0 = ps[ips];
    double force[3] = { 0,0,0 };
    for(unsigned int jps=0;jps<ps.size();jps++){
      if ( ips == jps ) continue;
      const SParticle& ps1 = ps[jps];
      Force(force,
          ps0,ps1,sph_radcutoff,SpikyKern,LapKern * SPH_VISC);
    }
    ps0.f[0] = force[0];
    ps0.f[1] = force[1];
    ps0.f[2] = force[2];
  }
}

template <class Search>
void SPH_AccumulateForceHash(
    std::vector<SParticle>& ps,
    double sph_radcutoff,
    double SPH_VISC,
    const Search& sg)
{
//  const double Poly6Kern = 315.0 / ( 64.0 * M_PI * pow( H, 9 ) );
  const double SpikyKern = -45.0 / ( M_PI * pow( sph_radcutoff, 6 ) );
  const double LapKern   = 45.0 / ( M_PI * pow( sph_radcutoff, 6 ) );
  std::vector<unsigned int> aIP;
  for(unsigned int ips=0;ips<ps.size();ips++){
    SParticle& ps0 = ps[ips];
    double force[3] = { 0,0,0 };
    sg.GetOneRingNeighbor(aIP,ps0.r.p);
    for(unsigned int jps : aIP){
      if ( ips == jps ) continue;
      const SParticle& ps1 = ps[jps];
      Force(force,
          ps0,ps1,sph_radcutoff,SpikyKern,LapKern * SPH_VISC);
    }
    ps0.f[0] = force[0];
    ps0.f[1] = force[1];
    ps0.f[2] = force[2];
  }
}

void SPH_UpdatePosition(
    std::vector<SParticle>& ps,
    double SPH_PMASS,
    double SPH_LIMIT,
    double SPH_EXTSTIFF,
    double SPH_EXTDAMP,
    double DT,
    double EPSILON,
    const double bbmin[3],
    const double bbmax[3],
    double radcutoff)
{
  const double g[3] = { 0.0, -9.8, 0.0};
  for(auto & p : ps){
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
        const double diff = 2.0 * radcutoff - (p.r[idim] - bbmin[idim] );
        if ( diff > EPSILON ){
          double norm[3] = { 0.0, 0.0, 0.0 }; norm[idim] = +1.0;
          const double adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * dfm2::Dot3( norm, p.v.p );
          accel[0] += adj * norm[0];
          accel[1] += adj * norm[1];
          accel[2] += adj * norm[2];
        }
      }
      {
        const double diff = 2.0 * radcutoff - ( bbmax[idim] - p.r[idim] );
        if ( diff > EPSILON ){
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
    p.r[0] += p.v[0] * DT;
    p.r[1] += p.v[1] * DT;
    p.r[2] += p.v[2] * DT;
  }
}

// --------------------------------

void myGlutDisplay(
    const std::vector<SParticle>& ps,
    const double MIN[3],
    const double MAX[3])
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
  for(auto & p : ps){
    ::glVertex3dv(p.r.p);
  }
  ::glEnd();
  delfem2::opengl::DrawBox3_Edge(MIN,MAX);
}


int main()
{
  const double DT = 0.001;
  const double sph_radcutoff = 0.01; // radius of influence
  const double sph_radius = 0.001; // radius of particle for collision to the walls
  const double EPSILON = 1.0e-5;
  const double bbmin[3] = { 0.0, 0.0,  0.0};
  const double bbmax[3] = { 0.1, 0.1,  0.05};
//  const double SPH_PMASS = 0.00020543; // particle mass (kg)
  const double SPH_PMASS = 0.0001; // particle mass (kg)
  const double SPH_INTSTIFF = 1.00;  //
  const double SPH_EXTSTIFF = 50000.0; // penalty coefficient of the wall repulsion force
  const double SPH_EXTDAMP  = 256.0; // damping of the wall repulsion force
  const double SPH_RESTDENSITY = 600.0; // kg / m^3
  const double SPH_VISC = 0.2; // pascal-second (Pa.s) = 1 kg m^-1 s^-1
  const double SPH_LIMIT = 200.0; // maximum speed

  std::vector<SParticle> ps;
  {
    std::mt19937 rndeng( std::random_device{}() );
//    std::mt19937 rndeng( 1 );
    std::uniform_real_distribution<double> distm0p1(-1., +1.);
    double d = pow( SPH_PMASS / SPH_RESTDENSITY, 1/3.0 ) * 0.87;
    const double INITMIN[3] = { 0.0, 0.0, 0.0};
    const double INITMAX[3] = { 0.05, 0.1, 0.05};
    const int ndivx = static_cast<int>(ceil((INITMAX[0]-INITMIN[0])/d));
    const int ndivy = static_cast<int>(ceil((INITMAX[1]-INITMIN[1])/d));
    const int ndivz = static_cast<int>(ceil((INITMAX[2]-INITMIN[2])/d));
    for (int idivx=0;idivx<ndivx;++idivx) {
      for (int idivy=0;idivy<ndivy;++idivy) {
        for (int idivz=0;idivz<ndivz;++idivz) {
          SParticle p;
          p.r[0] = INITMIN[0] + d*(idivx+0.5+distm0p1(rndeng)*0.3);
          p.r[1] = INITMIN[1] + d*(idivy+0.5+distm0p1(rndeng)*0.3);
          p.r[2] = INITMIN[2] + d*(idivz+0.5+distm0p1(rndeng)*0.3);
          // rejection sampling
          if( p.r[0] < INITMIN[0] || p.r[0] > INITMAX[0] ){ continue; }
          if( p.r[1] < INITMIN[1] || p.r[1] > INITMAX[1] ){ continue; }
          if( p.r[2] < INITMIN[2] || p.r[2] > INITMAX[2] ){ continue; }
          p.v[0] = 0.0;
          p.v[1] = 0.0;
          p.v[2] = 0.0;
          ps.push_back( p );
        }
      }
    }
  }
  std::cout << "particle size : " << ps.size() << std::endl;

  dfm2::glfw::CViewer3 viewer(0.1);
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  {
    ::glfwSetWindowTitle(viewer.window, "SPH with Spatial Hash");
    delfem2::SearchGrid sg;
    sg.Initialize(bbmin, bbmax, sph_radcutoff, ps.size());
    for (unsigned int ip = 0; ip < ps.size(); ++ip) {
      sg.aGrid2Obj[ip].igrid = sg.GetGridIndex(ps[ip].r.p);
      sg.aGrid2Obj[ip].iobj = ip;
    }
    sg.PostProcess(true);
    for (unsigned int iframe = 0; iframe < 300; ++iframe) {
      for(auto& gridobj : sg.aGrid2Obj){
        const unsigned int ip = gridobj.iobj;
        gridobj.igrid = sg.GetGridIndex(ps[ip].r.p);
      }
      sg.PostProcess(false);
      SPH_DensityPressureHash(
          ps,
          sph_radcutoff, SPH_PMASS, SPH_RESTDENSITY, SPH_INTSTIFF, sg);
      SPH_AccumulateForceHash(
          ps,
          sph_radcutoff, SPH_VISC, sg);
      SPH_UpdatePosition(
          ps,
          SPH_PMASS, SPH_LIMIT, SPH_EXTSTIFF, SPH_EXTDAMP,
          DT, EPSILON, bbmin, bbmax, sph_radius);
      // -----
      viewer.DrawBegin_oldGL();
      myGlutDisplay(ps, bbmin, bbmax);
      viewer.SwapBuffers();
      glfwPollEvents();
      if (glfwWindowShouldClose(viewer.window)) {
        glfwDestroyWindow(viewer.window);
        glfwTerminate();
        exit(EXIT_SUCCESS);
      }
    }
  }

  // --------------------
  ::glfwSetWindowTitle(viewer.window, "SPH without Spatial Hash");
  for(unsigned int iframe=0;iframe<30;++iframe){
    SPH_DensityPressure(
        ps,
        sph_radcutoff, SPH_PMASS, SPH_RESTDENSITY, SPH_INTSTIFF);
    SPH_AccumulateForce(
        ps,
        sph_radcutoff, SPH_VISC);
    SPH_UpdatePosition(
        ps,
        SPH_PMASS,SPH_LIMIT,SPH_EXTSTIFF,SPH_EXTDAMP,
        DT,EPSILON,bbmin,bbmax,sph_radius);
    // -----
    viewer.DrawBegin_oldGL();
    myGlutDisplay(ps,bbmin,bbmax);
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){
      glfwDestroyWindow(viewer.window);
      glfwTerminate();
      exit(EXIT_SUCCESS);
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
