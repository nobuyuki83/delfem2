/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <iostream>
#include <vector>
#include <cmath>
#include <random>
// ------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"

template <typename VAL>
VAL max(VAL i, VAL j) { return (i>j?i:j); }

template <typename VAL>
VAL min(VAL i, VAL j) { return (i>j?j:i); }

// Clamped Fetch
double ClampedFetch(
    double* x,
    unsigned int i, 
    unsigned int j,
    unsigned int nw,
    unsigned int nh)
{
	i = min(max((unsigned int)0,i),nw-1);
	j = min(max((unsigned int)0,j),nh-1);
	return x[i+j*nw];
}

// Gauss-Seidel Iteration
void GaussSeidelSolver(
    double* p,
    const double* d,
    unsigned int n,
    unsigned int K0 )
{
	double h2 = 1.0/(n*n);
	for(unsigned int k=0; k<K0; k++ ){
    double t = 0;
		for(unsigned int i=0; i<n; i++ ){
      for(unsigned int j=0; j<n; j++ ){
        const double p0 = p[i+j*n];
        p[i+j*n] = (
            ClampedFetch(p,i+1,j,n,n) +
            ClampedFetch(p,i-1,j,n,n) +
            ClampedFetch(p,i,j+1,n,n) +
            ClampedFetch(p,i,j-1,n,n) - h2*d[i+j*n] ) / 4.0;
        const double p1 = p[i+j*n];
        t += fabs(p0-p1);
      }
    }
    if( k % 100 == 0 ){
      std::cout << "itr : " << k << " " << t << std::endl;
    }
	}
}

void CompDivergence(
    std::vector<double>& divag,
    unsigned int ngrid,
    const std::vector<double>& velou,
    const std::vector<double>& velov)
{
  for(unsigned int ig=0;ig<ngrid;ig++){
    for(unsigned int jg=0;jg<ngrid;jg++){
      double div = (velou[ig+1+jg*(ngrid+1)]-velou[ig+jg*(ngrid+1)])
                 + (velov[ig+(jg+1)*ngrid]-velov[ig+jg*ngrid]);
      divag[ig+jg*ngrid] = div*ngrid;
    }
	}
}

void EnforceBoundary(
    std::vector<double>& velou,
    std::vector<double>& velov,
    unsigned int ngrid)
{
  for(unsigned int ig=0;ig<ngrid;ig++){ velov[ig] = 0; }
  for(unsigned int ig=0;ig<ngrid;ig++){ velov[ig+ngrid*ngrid] = 0; }  
  for(unsigned int jg=0;jg<ngrid;jg++){ 
    velou[      jg*(ngrid+1)] = 0;
    velou[ngrid+jg*(ngrid+1)] = 0;    
  }
}

void CompPressure(
    std::vector<double>& divag,
    std::vector<double>& press,
    unsigned int ngrid,
    double rho,
    double dt)
{
  const double dtmp1 = rho/dt;
  for(unsigned int ig=0;ig<ngrid;ig++){
    for(unsigned int jg=0;jg<ngrid;jg++){
      divag[ig+jg*ngrid] *= dtmp1;
    }
  }
  GaussSeidelSolver(press.data(), divag.data(), ngrid, 1000);
}

void SubtractPressure(
    std::vector<double>& velou,
    std::vector<double>& velov,
    unsigned int ngrid,
    double dt,
    double rho,
    const std::vector<double>& press)
{
  double dtmp1 = dt/rho;
  for(unsigned int ig=0;ig<ngrid-1;ig++){
    for(unsigned int jg=0;jg<ngrid;jg++){
      velou[(ig+1)+jg*(ngrid+1)] -= dtmp1*(press[(ig+1)+jg*ngrid]-press[ig+jg*ngrid])*ngrid;
    }
  }
  for(unsigned int ig=0;ig<ngrid;ig++){
    for(unsigned int jg=0;jg<ngrid-1;jg++){
      velov[ig+(jg+1)*ngrid] -= dtmp1*(press[ig+(jg+1)*ngrid]-press[ig+(jg+0)*ngrid])*ngrid;
    }
  }  
}


double linear_interpolate (
    const double* d,
    unsigned int ndim,
    unsigned int idim,
    unsigned int nw,
    unsigned int nh,
    double x,
    double y )
{
	x = max(0.0,min((double)nw,x));
	y = max(0.0,min((double)nh,y));
	int i = min(x,(double)nw-2);
	int j = min(y,(double)nh-2);
	return ( (i+1-x)*d[(i+j*nw    )*ndim+idim]+(x-i)*d[(i+1+j*nw    )*ndim+idim])*(j+1-y) 
        + ((i+1-x)*d[(i+(j+1)*nw)*ndim+idim]+(x-i)*d[(i+1+(j+1)*nw)*ndim+idim])*(y-j);
}

void CompAdvectionSemiLagrangian(
    std::vector<double>& velou,
    std::vector<double>& velov,
    std::vector<double>& velou_tmp,
    std::vector<double>& velov_tmp,
    unsigned int ngrid,
    double dt)
{  
  for(unsigned int jg=0;jg<ngrid+0;jg++){
    for(unsigned int ig=0;ig<ngrid+1;ig++){
      velou_tmp[(ig+jg*(ngrid+1))*2+0] = velou[ig+jg*(ngrid+1)];
      velou_tmp[(ig+jg*(ngrid+1))*2+1] = ( ClampedFetch(velov.data(), ig-1,jg,   ngrid,ngrid+1)
                                          +ClampedFetch(velov.data(), ig,  jg,   ngrid,ngrid+1)
                                          +ClampedFetch(velov.data(), ig-1,jg+1, ngrid,ngrid+1)
                                          +ClampedFetch(velov.data(), ig,  jg+1, ngrid,ngrid+1) )*0.25;
    }
	}
  
  for(unsigned int jg=0;jg<ngrid+1;jg++){
    for(unsigned int ig=0;ig<ngrid+0;ig++){
      velov_tmp[(ig+jg*ngrid)*2+0] = ( ClampedFetch(velou.data(), ig,  jg-1,ngrid+1,ngrid)
                                      +ClampedFetch(velou.data(), ig+1,jg-1,ngrid+1,ngrid)
                                      +ClampedFetch(velou.data(), ig,  jg,  ngrid+1,ngrid)
                                      +ClampedFetch(velou.data(), ig+1,jg,  ngrid+1,ngrid) )*0.25;
      velov_tmp[(ig+jg*ngrid)*2+1] = velov[ig+jg*ngrid];
    }
	}

  for(unsigned int jg=0;jg<ngrid+0;jg++){
    for(unsigned int ig=0;ig<ngrid+1;ig++){
      const double* velo = velou_tmp.data()+(ig+jg*(ngrid+1))*2;
      const double p[2] = {ig-velo[0]*dt*ngrid,jg-velo[1]*dt*ngrid};
      const double u = linear_interpolate(velou_tmp.data(),2,0,ngrid+1,ngrid,p[0],p[1]);
      velou[ig+jg*(ngrid+1)] = u;
    }
  }

  for(unsigned int jg=0;jg<ngrid+1;jg++){
    for(unsigned int ig=0;ig<ngrid+0;ig++){
      const double* velo = velov_tmp.data()+(ig+jg*ngrid)*2;
      const double p[2] = {ig-velo[0]*dt*ngrid,jg-velo[1]*dt*ngrid};
      const double v = linear_interpolate(velov_tmp.data(),2,1,ngrid,ngrid+1,p[0],p[1]);
      velov[ig+jg*ngrid] = v;
    }
  }
}


void AssignGravity(
    std::vector<double>& velou,
    std::vector<double>& velov,
    unsigned int ngrid,
    const double gravity[2],
    double dt)
{
  for(unsigned int jg=0;jg<ngrid+0;jg++){
    for(unsigned int ig=0;ig<ngrid+1;ig++){
      velou[ig+jg*(ngrid+1)] += gravity[0]*dt;
    }
  }
  
  for(unsigned int jg=0;jg<ngrid+1;jg++){
    for(unsigned int ig=0;ig<ngrid+0;ig++){        
      velov[ig+jg*ngrid] += gravity[1]*dt;
    }
  }  
  
}

void glutMyDisplay(
    unsigned int ngrid,
    const std::vector<double>& velou,
    const std::vector<double>& velov,
    const std::vector<double>& press)
{
	glClear(GL_COLOR_BUFFER_BIT);
  {
    double h = 1.0/ngrid;    
    ::glBegin(GL_QUADS);
    for(unsigned int ig=0;ig<ngrid;ig++){
      for(unsigned int jg=0;jg<ngrid;jg++){
        double p = press[ig+jg*ngrid];
        //			glColor4d(p>0,0.0,p<0,fabs(p));    
        glColor4d(p>0,0.0,p<0,0.8);
        ::glVertex2d((ig+0)*h,(jg+0)*h);
        ::glVertex2d((ig+1)*h,(jg+0)*h);
        ::glVertex2d((ig+1)*h,(jg+1)*h);      
        ::glVertex2d((ig+0)*h,(jg+1)*h);
      }
    }
    ::glEnd();    
  }

  { // draw velocity
    double h = 1.0/ngrid;
    ::glColor3d(0,1,1);
    ::glBegin(GL_LINES);
    for(unsigned int ig=0;ig<ngrid;ig++){
      for(unsigned int jg=0;jg<ngrid;jg++){
        const double p[2] = {(ig+0.5)*h,(jg+0.5)*h};
        const double u = (velou[ig+jg*(ngrid+1)] + velou[ig+1+ jg   *(ngrid+1)])*0.5*3;
        const double v = (velov[ig+jg*ngrid    ] + velov[ig+  (jg+1)* ngrid   ])*0.5*3;
        ::glVertex2d(p[0],  p[1]  );
        ::glVertex2d(p[0]+u,p[1]+v);
      }
    }
    ::glEnd();
  }
}


int main (int argc, char * argv[])
{
  const unsigned int ngrid = 32;
  std::vector<double> velou(ngrid*(ngrid+1), 0.0);  // (ngrid+1)*ngrid
  std::vector<double> velov(ngrid*(ngrid+1), 0.0);  // ngrid*(ngrid+1)
  std::vector<double> press(ngrid*ngrid, 0.0);  // ngrid*ngrid
  std::vector<double> divag(ngrid*ngrid, 0.0);  // ngrid*ngrid;
  std::vector<double> velou_tmp;  // (ngrid+1)*ngrid*2
  std::vector<double> velov_tmp;  // ngrid*(ngrid+1)*2
  // -----------
  // allocate memory
  velou_tmp.resize(ngrid*(ngrid+1)*2);
  velov_tmp.resize(ngrid*(ngrid+1)*2);
  // -----------
  const double dt = 0.1;
  const double rho = 1.0;
  const double gravity[2] = {0,-0.01};
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  
  int iframe = -1;
  while (true)
  {
    iframe = (iframe+1)%100;
    if( iframe == 0 ){
      for(unsigned int i=0;i<ngrid*(ngrid+1);i++){ velou[i] = dist(mt); }
      for(unsigned int i=0;i<ngrid*(ngrid+1);i++){ velov[i] = dist(mt); }
    }
    // ----
    AssignGravity(velou,velov,ngrid,gravity,dt);
    EnforceBoundary(velou,velov, ngrid);
    CompDivergence(divag, ngrid, velou, velov);
    CompPressure(divag,press,
        ngrid,rho,dt);
    SubtractPressure(velou,velov,
        ngrid,dt,rho,press);
    CompAdvectionSemiLagrangian(velou,velov,velou_tmp,velov_tmp,
        ngrid,dt);
    // ----
    viewer.DrawBegin_oldGL();
    glutMyDisplay(ngrid,velou,velov,press);
    viewer.DrawEnd_oldGL();
    if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
