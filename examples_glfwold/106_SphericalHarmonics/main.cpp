/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include "delfem2/vec3.h"
#include "delfem2/specialfuncs.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -----------------------------------

static void drawShphere_Heatmap
(double (*value)(double,double,double),
 void (*color)(double) )
{
  const double pi = 3.1415926535;
  int nla = 32; double dl = pi/nla;
  int nlo = 64; double dr = 2.0*pi/nlo;
  ::glBegin(GL_QUADS);
  for(int ila=0;ila<nla-1;ila++){
    int ila0 = ila;
    int ila1 = ila+1;
    double y0 = cos(dl*ila0);
    double y1 = cos(dl*ila1);
    double r0 = sin(dl*ila0);
    double r1 = sin(dl*ila1);
    for(int ilo=0;ilo<nlo;ilo++){
      int ilo0 = ilo;
      int ilo1 = ilo+1;
      double x0 = sin(dr*ilo0);
      double x1 = sin(dr*ilo1);
      double z0 = cos(dr*ilo0);
      double z1 = cos(dr*ilo1);
      dfm2::CVec3d a(r0*x0,y0,r0*z0);
      dfm2::CVec3d b(r0*x1,y0,r0*z1);
      dfm2::CVec3d c(r1*x1,y1,r1*z1);
      dfm2::CVec3d d(r1*x0,y1,r1*z0);
      {
        double x = (r0+r1)*(x0+x1)*0.25;
        double y = (y0+y1)*0.5;
        double z = (r0+r1)*(z0+z1)*0.25;
        double invl = 1.0/sqrt(x*x+y*y+z*z);
        x *= invl;
        y *= invl;
        z *= invl;
        glNormal3d(x,y,z);
      }
      double va = value(a.x(),a.y(),a.z());
      double vb = value(b.x(),b.y(),b.z());
      double vc = value(c.x(),c.y(),c.z());
      double vd = value(d.x(),d.y(),d.z());
      color(va+0.5); delfem2::opengl::myGlVertex(a);
      color(vb+0.5); delfem2::opengl::myGlVertex(b);
      color(vc+0.5); delfem2::opengl::myGlVertex(c);
      color(vd+0.5); delfem2::opengl::myGlVertex(d);
    }
  }
  ::glEnd();
}

static void drawShphere_Radius
(double (*value)(double,double,double),
 void (*color)(double) )
{
  const double pi = 3.1415926535;
  int nla = 64; double dla = pi/nla;
  int nlo =128; double dlo = 2.0*pi/nlo;
  ::glBegin(GL_TRIANGLES);
  for(int ila=0;ila<nla-1;ila++){
    int ila0 = ila;
    int ila1 = ila+1;
    double y0 = cos(dla*ila0);
    double y1 = cos(dla*ila1);
    double r0 = sin(dla*ila0);
    double r1 = sin(dla*ila1);
    for(int ilo=0;ilo<nlo;ilo++){
      int ilo0 = ilo;
      int ilo1 = ilo+1;
      double x0 = sin(dlo*ilo0);
      double x1 = sin(dlo*ilo1);
      double z0 = cos(dlo*ilo0);
      double z1 = cos(dlo*ilo1);
      dfm2::CVec3d a(r0*x0,y0,r0*z0);
      dfm2::CVec3d b(r0*x1,y0,r0*z1);
      dfm2::CVec3d c(r1*x1,y1,r1*z1);
      dfm2::CVec3d d(r1*x0,y1,r1*z0);
      double va = value(a.x(),a.y(),a.z());
      double vb = value(b.x(),b.y(),b.z());
      double vc = value(c.x(),c.y(),c.z());
      double vd = value(d.x(),d.y(),d.z());
      delfem2::opengl::myGlNormal(va*a,vc*c,vb*b);
      color(va+0.5); delfem2::opengl::myGlVertex(fabs(va)*a);
      color(vb+0.5); delfem2::opengl::myGlVertex(fabs(vb)*b);
      color(vc+0.5); delfem2::opengl::myGlVertex(fabs(vc)*c);
      delfem2::opengl::myGlNormal(vd*d,vc*c,va*a);
      color(vd+0.5); delfem2::opengl::myGlVertex(fabs(vd)*d);
      color(vc+0.5); delfem2::opengl::myGlVertex(fabs(vc)*c);
      color(va+0.5); delfem2::opengl::myGlVertex(fabs(va)*a);
    }
  }
  ::glEnd();
}


// -------------------------------

const int nl = 10;
int l=0;
int m=0;
std::vector<std::complex<double>> aDecmpSH;
std::vector<double> aDecmpSHReal;

// -------------------------------

double evaluateSH(double x, double y, double z)
{
  double a[100];
//  makeEvaluationArraySH(l+1,a, x,y,z);
  dfm2::makeArray_SphericalHarmonics(a, l+1, x,y,z);
  int ish = l*(l+1)+m;
  return a[ish];
}

int main(int argc,char* argv[])
{
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  
  // -----------------------------
  
  viewer.nav.camera.view_height = 1.0;
  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      ++m;
      if( m > l ){ l++; m=-l; }
      if( l >= nl ){ l=0; m=0; }
    }
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    drawShphere_Radius(evaluateSH,delfem2::opengl::heatmap_glDiffuse);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


