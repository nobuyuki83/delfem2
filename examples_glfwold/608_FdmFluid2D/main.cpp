#include <iostream>
#include <vector>
#include <math.h>
// ------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"

#define max(i,j) (i>j?i:j)
#define min(i,j) (i>j?j:i)

// Clamped Fetch
double ClampedFetch(
    double* x,
    int i, int j, unsigned int nw, unsigned int nh)
{
	i = min(max(0,i),nw-1);
	j = min(max(0,j),nh-1);
	return x[i+j*nw];
}

// Gauss-Seidel Iteration
void GaussSeidelSolver(
    double* p,
    const double* d,
    int n, int t )
{
	double h2 = 1.0/(n*n);
	for(unsigned int k=0; k<t; k++ ){
    double t = 0;
		for(int i=0; i<n; i++ ){
    for(int j=0; j<n; j++ ){
      const double p0 = p[i+j*n];
			p[i+j*n] = (ClampedFetch(p,i+1,j,n,n)+ClampedFetch(p,i-1,j,n,n)+ClampedFetch(p,i,j+1,n,n)+ClampedFetch(p,i,j-1,n,n)
                  -h2*d[i+j*n]) / 4.0;
      const double p1 = p[i+j*n];
      t += fabs(p0-p1);
		}
    }
    if( k % 100 == 0 ){
      std::cout << "itr : " << k << " " << t << std::endl;
    }
	}
}


const unsigned int ngrid = 32;
std::vector<double> velou;  // (ngrid+1)*ngrid
std::vector<double> velov;  // ngrid*(ngrid+1)
std::vector<double> press;  // ngrid*ngrid
std::vector<double> divag;  // ngrid*ngrid;
const double dt = 0.1;
const double rho = 1;
bool is_animation = true;
const double gravity[2] = {0,-0.01};

std::vector<double> velou_tmp;  // (ngrid+1)*ngrid*2
std::vector<double> velov_tmp;  // ngrid*(ngrid+1)*2

std::vector<double> aMarker;

void ClearValue(){
  for(unsigned int i=0;i<ngrid*(ngrid+1);i++){ velou[i] = 0; }
  for(unsigned int i=0;i<ngrid*(ngrid+1);i++){ velov[i] = 0; }  
  for(unsigned int i=0;i<ngrid*ngrid;i++){ press[i] = 0; }
  for(unsigned int i=0;i<ngrid*ngrid;i++){ divag[i] = 0; }  
  aMarker.clear();
}

void InitSimulation()
{
  velou.resize(ngrid*(ngrid+1));
  velov.resize(ngrid*(ngrid+1));
  press.resize(ngrid*ngrid);
  divag.resize(ngrid*ngrid);
  velou_tmp.resize(ngrid*(ngrid+1)*2);
  velov_tmp.resize(ngrid*(ngrid+1)*2);
  ClearValue();
}

void CompDivergence() {
  for(unsigned int ig=0;ig<ngrid;ig++){
  for(unsigned int jg=0;jg<ngrid;jg++){    
		double div = (velou[ig+1+jg*(ngrid+1)]-velou[ig+jg*(ngrid+1)]) 
               + (velov[ig+(jg+1)*ngrid]-velov[ig+jg*ngrid]);
		divag[ig+jg*ngrid] = div*ngrid;
  }
	}
}

void EnforceBoundary(){
  for(unsigned int ig=0;ig<ngrid;ig++){ velov[ig] = 0; }
  for(unsigned int ig=0;ig<ngrid;ig++){ velov[ig+ngrid*ngrid] = 0; }  
  for(unsigned int jg=0;jg<ngrid;jg++){ 
    velou[      jg*(ngrid+1)] = 0;
    velou[ngrid+jg*(ngrid+1)] = 0;    
  }
}

void CompPressure(){
  const double dtmp1 = rho/dt;
  for(unsigned int ig=0;ig<ngrid;ig++){
  for(unsigned int jg=0;jg<ngrid;jg++){     
		divag[ig+jg*ngrid] *= dtmp1;    
  }
  }
  GaussSeidelSolver(press.data(), divag.data(), ngrid, 1000);
}

void SubtractPressure(){
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


static double linear_interpolate (double* d, unsigned int ndim, unsigned int idim, int nw, int nh, double x, double y ) 
{
	x = max(0.0,min(nw,x));
	y = max(0.0,min(nh,y));
	int i = min(x,nw-2);
	int j = min(y,nh-2);
	return ( (i+1-x)*d[(i+j*nw    )*ndim+idim]+(x-i)*d[(i+1+j*nw    )*ndim+idim])*(j+1-y) 
        + ((i+1-x)*d[(i+(j+1)*nw)*ndim+idim]+(x-i)*d[(i+1+(j+1)*nw)*ndim+idim])*(y-j);
}

void CompAdvectionSemiLagrangian()
{  
  for(int jg=0;jg<ngrid+0;jg++){
  for(int ig=0;ig<ngrid+1;ig++){    
		velou_tmp[(ig+jg*(ngrid+1))*2+0] = velou[ig+jg*(ngrid+1)];
		velou_tmp[(ig+jg*(ngrid+1))*2+1] = ( ClampedFetch(velov.data(), ig-1,jg,   ngrid,ngrid+1)
                                        +ClampedFetch(velov.data(), ig,  jg,   ngrid,ngrid+1)
                                        +ClampedFetch(velov.data(), ig-1,jg+1, ngrid,ngrid+1)
                                        +ClampedFetch(velov.data(), ig,  jg+1, ngrid,ngrid+1) )*0.25;
  }
	}
  
  for(int jg=0;jg<ngrid+1;jg++){    
  for(int ig=0;ig<ngrid+0;ig++){    
    velov_tmp[(ig+jg*ngrid)*2+0] = ( ClampedFetch(velou.data(), ig,  jg-1,ngrid+1,ngrid)
                                    +ClampedFetch(velou.data(), ig+1,jg-1,ngrid+1,ngrid)
                                    +ClampedFetch(velou.data(), ig,  jg,  ngrid+1,ngrid)
                                    +ClampedFetch(velou.data(), ig+1,jg,  ngrid+1,ngrid) )*0.25;
    velov_tmp[(ig+jg*ngrid)*2+1] = velov[ig+jg*ngrid];    
  }
	}
  
  {
    const unsigned int nMark = aMarker.size()/2;
    for(unsigned int imark=0;imark<nMark;imark++){
      const double q0[2] = {aMarker[imark*2+0]*ngrid,aMarker[imark*2+1]*ngrid};
      const double u = linear_interpolate(velou_tmp.data(),2,0,ngrid+1,ngrid,q0[0],q0[1]);
      const double v = linear_interpolate(velov_tmp.data(),2,1,ngrid,ngrid+1,q0[0],q0[1]);
      const double p0[2] = {aMarker[imark*2+0],aMarker[imark*2+1]};
      double p1[2] = {p0[0]+u*dt,p0[1]+v*dt};
      if( p1[0] < 0 ){ p1[0] = -p1[0]; }
      if( p1[0] > 1 ){ p1[0] = 2-p1[0]; }      
      if( p1[1] < 0 ){ p1[1] = -p1[1]; }
      if( p1[1] > 1 ){ p1[1] = 2-p1[1]; }            
      aMarker[imark*2+0] = p1[0];
      aMarker[imark*2+1] = p1[1];
    }
  }
  
  for(int jg=0;jg<ngrid+0;jg++){
  for(int ig=0;ig<ngrid+1;ig++){    
    const double* velo = velou_tmp.data()+(ig+jg*(ngrid+1))*2;
    const double p[2] = {ig-velo[0]*dt*ngrid,jg-velo[1]*dt*ngrid};    
    const double u = linear_interpolate(velou_tmp.data(),2,0,ngrid+1,ngrid,p[0],p[1]);
    velou[ig+jg*(ngrid+1)] = u;
  }
  }

  for(int jg=0;jg<ngrid+1;jg++){
  for(int ig=0;ig<ngrid+0;ig++){    
    const double* velo = velov_tmp.data()+(ig+jg*ngrid)*2;
    const double p[2] = {ig-velo[0]*dt*ngrid,jg-velo[1]*dt*ngrid};    
    const double v = linear_interpolate(velov_tmp.data(),2,1,ngrid,ngrid+1,p[0],p[1]);
    velov[ig+jg*ngrid] = v;
  }
  }
}


void AssignGravity(){
  for(int jg=0;jg<ngrid+0;jg++){
    for(int ig=0;ig<ngrid+1;ig++){    
      velou[ig+jg*(ngrid+1)] += gravity[0]*dt;
    }
  }
  
  for(int jg=0;jg<ngrid+1;jg++){
    for(int ig=0;ig<ngrid+0;ig++){        
      velov[ig+jg*ngrid] += gravity[1]*dt;
    }
  }  
  
}

  
void StepTime(){
  AssignGravity();
  EnforceBoundary();
  CompDivergence();
  CompPressure();
  SubtractPressure();
  CompAdvectionSemiLagrangian();
}

void glutMyDisplay(void){
  if( is_animation ){
    StepTime();
  }
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
  /*
  {
    double h = 1.0/ngrid;
    ::glBegin(GL_POINTS);
    for(unsigned int ig=0;ig<ngrid;ig++){
    for(unsigned int jg=0;jg<ngrid;jg++){
      const double p[2] = {(ig+0.5)*h,(jg+0.5)*h};    
      ::glVertex2d(p[0],p[1]);
    }    
    }
    ::glEnd();
  }
   */
  {
    const unsigned int nMark = aMarker.size()/2;
    ::glColor3d(1,0,0);
    ::glPointSize(5);
    ::glBegin(GL_POINTS);
    for(unsigned int imark=0;imark<nMark;imark++){
      const double p[2] = {aMarker[imark*2+0],aMarker[imark*2+1]};           
      ::glVertex3d(p[0],p[1],0.1);
    }    
    ::glEnd();
  }
}

/*
static void init( int gsize ) {
	glClearColor(0.0, 0.0, 0.0, 1.0);
}
 */

/*
static void glutMyReshape(int w, int h) {
  ::glViewport(0,0,w,h);
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::glOrtho(0,1, 0,1, -1,1);
}

static void glutMyIdle( void ) {

	glutPostRedisplay();
}

static void glutMyKeyboard( unsigned char key, int x, int y ) {
  if( key == ' ' ){
    ClearValue();
  }
  if( key == 'b' ){
    const double h = 1.0/ngrid;
    for(unsigned int jg=0;jg<ngrid;jg++){
    for(unsigned int ig=0;ig<ngrid+1;ig++){    
      const double p[2] = {ig*h,(jg+0.5)*h};
      const double r = sqrt( (p[0]-0.5)*(p[0]-0.5) + (p[1]-0.5)*(p[1]-0.5) );
      velou[ig+jg*(ngrid+1)] += max(0,0.001-r*r);
    }
    }
  }
  if( key == 'c' ){
    StepTime();
  }
  else if( key == 'a' ){
    is_animation = !is_animation;
  }
}

int prev_x, prev_y;
int imodifier;

static void glutMyMouse ( int button, int state, int x, int y ) {
	prev_x = x;
	prev_y = y;
  imodifier = glutGetModifiers();
  if( imodifier == GLUT_ACTIVE_SHIFT ){
    int w[4];
    ::glGetIntegerv(GL_VIEWPORT,w);
    int win_x = w[2];
    int win_y = w[3];        
    const double px = x/(GLdouble)win_x;
    const double py = 1.0 - y/(GLdouble)win_y;
    aMarker.push_back(px);
    aMarker.push_back(py);
  }
}
 */

void SetVeloFieldAround(double x, double y, double u, double v){
  const double h = 1.0/ngrid;
  for(unsigned int jg=0;jg<ngrid;jg++){
  for(unsigned int ig=0;ig<ngrid+1;ig++){    
    const double p[2] = {ig*h,(jg+0.5)*h};
    const double r = sqrt( (p[0]-x)*(p[0]-x) + (p[1]-y)*(p[1]-y) );
    velou[ig+jg*(ngrid+1)] += max(0,0.01-r*r)*u*100;
  }
  }  
  for(unsigned int jg=0;jg<ngrid+1;jg++){
  for(unsigned int ig=0;ig<ngrid;ig++){    
    const double p[2] = {(ig+0.5)*h,jg*h};
    const double r = sqrt( (p[0]-x)*(p[0]-x) + (p[1]-y)*(p[1]-y) );
    velov[ig+jg*ngrid] += max(0,0.01-r*r)*v*100;
  }
  }    
}

/*
static void glutMyMotion ( int x, int y ) {
  if( imodifier == GLUT_ACTIVE_SHIFT ){
    int w[4];
    ::glGetIntegerv(GL_VIEWPORT,w);
    int win_x = w[2];
    int win_y = w[3];        
    const double px = x/(GLdouble)win_x;
    const double py = 1.0 - y/(GLdouble)win_y;
    aMarker.push_back(px);
    aMarker.push_back(py);
  }
  else{
    int w[4];
    ::glGetIntegerv(GL_VIEWPORT,w);
    int win_x = w[2];
    int win_y = w[3];    
    SetVeloFieldAround(x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y,
                       (x-prev_x)/(GLdouble)win_x, -(y-prev_y)/(GLdouble)win_y );
	}
	prev_x = x;
	prev_y = y;
}
 */


int main (int argc, char * argv[])
{
  InitSimulation();
  ClearValue();
  
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  
  int iframe = -1;
  while (true)
  {
    iframe = (iframe+1)%100;
    if( iframe == 0 ){
      for(unsigned int i=0;i<ngrid*(ngrid+1);i++){ velou[i] = 2.0*rand()/(RAND_MAX+1.0)-1.0; }
      for(unsigned int i=0;i<ngrid*(ngrid+1);i++){ velov[i] = 2.0*rand()/(RAND_MAX+1.0)-1.0; }
    }
    viewer.DrawBegin_oldGL();
    glutMyDisplay();
    viewer.DrawEnd_oldGL();
    if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
