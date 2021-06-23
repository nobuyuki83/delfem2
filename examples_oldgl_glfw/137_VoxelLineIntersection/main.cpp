/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/mshmisc.h"
#include "delfem2/gridvoxel.h"
#include <GLFW/glfw3.h>
#include <random>


namespace dfm2 = delfem2;

// ------------------------------------------------------

void Draw_CGrid3
(const dfm2::CGrid3<int>& grid,
 std::vector<unsigned int>& aIdVox)
{
  { // set-up transformation
    const dfm2::CMat4d& am = grid.am;
    dfm2::CMat4d amt = am.Transpose();
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    ::glMultMatrixd(amt.mat);
  }
  // -----------
  /*
  {
    ::glDisable(GL_LIGHTING);
    const double x0 = (double)grid.ndivx;
    const double y0 = (double)grid.ndivy;
    const double z0 = (double)grid.ndivz;
    double p[2][3] = {
      { 0, 0, 0},
      {x0, 0, 0},
    };
    ::glBegin(GL_LINES);
    ::glVertex3dv(p[0]);
    ::glVertex3dv(p[1]);
    ::glEnd();
  }
   */
  {
    dfm2::opengl::myGlMaterialDiffuse(dfm2::CColor::Gray(0.9));
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_NORMALIZE);
    const unsigned int nx = grid.ndivx;
    const unsigned int ny = grid.ndivy;
    const unsigned int nz = grid.ndivz;
    for(unsigned int iz=0;iz<nz;++iz){
      for(unsigned int iy=0;iy<ny;++iy){
        for(unsigned int ix=0;ix<nx;++ix){
          if( grid.aVal[iz*ny*nx+iy*nx+ix] == 0 ){ continue; }
          const double pmin[3] = {(double)ix,(double)iy,(double)iz};
          const double pmax[3] = {(double)ix+1.,(double)iy+1.,(double)iz+1.};
          dfm2::opengl::DrawBox3_Face(pmin,pmax);
        }
      }
    }
  }
  ::glClear(GL_DEPTH_BUFFER_BIT);
  
  {
    ::glDisable(GL_LIGHTING);
//    ::glDisable(GL_DEPTH_TEST);
    const unsigned int nx = grid.ndivx;
    const unsigned int ny = grid.ndivy;
    for(unsigned int iivox=0;iivox<aIdVox.size();++iivox){
      unsigned int ivox0 = aIdVox[iivox];
      const int iz0 = ivox0/(ny*nx);
      const int iy0 = (ivox0-iz0*ny*nx)/nx;
      const int ix0 = ivox0-iz0*ny*nx-iy0*nx;
      const double pmin[3] = {(double)ix0,(double)iy0,(double)iz0};
      const double pmax[3] = {(double)ix0+1.,(double)iy0+1.,(double)iz0+1.};
      ::glColor3d(1,0,0);
      dfm2::opengl::DrawBox3_Face(pmin, pmax);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawBox3_Edge(pmin, pmax);
//      dfm2::opengl::DrawSphereAt(32, 32, 1.0, ix0+0.5, iy0+0.5, iz0+0.5);
    }
  }
  
  ::glClear(GL_DEPTH_BUFFER_BIT);

  
  // --------
  // end transformation
  ::glPopMatrix();
}


// ------------------------------------------------------



int main(int argc,char* argv[])
{
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  viewer.camera.Rot_Camera(+0.2, -0.2);
  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);
  
  dfm2::CGrid3<int> grid;
  {
    const unsigned int n = 32;
    const double el = 1.0/n;
    grid.am.SetScale(el, el, el);
    grid.am.mat[ 3] = -0.5;
    grid.am.mat[ 7] = -0.5;
    grid.am.mat[11] = -0.5;
    grid.Initialize(n, n, n, 1);
  }
  std::random_device rd;
  std::mt19937 reng(rd());
  std::uniform_int_distribution<unsigned int> dist(0,grid.aVal.size()-1);
  
  // -------
  int iframe = 0;
  while (true)
  {
    std::vector<unsigned int> aIndVox;
    iframe++;
//    iframe = 1;
    dfm2::CVec3d ps(-0.3+0.9*sin(iframe*0.01),-0.3+0.9*sin(iframe*0.03),-0.3+0.9*sin(iframe*0.05));
    dfm2::CVec3d pe(+0.3+0.9*sin(iframe*0.07),+0.3+0.9*sin(iframe*0.11),+0.3+0.9*sin(iframe*0.13));
    Intersection_VoxelGrid_LinSeg(aIndVox,
                                  grid, ps, pe);
    // -------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_DEPTH_TEST);
    Draw_CGrid3(grid,aIndVox);
    ::glDisable(GL_DEPTH_TEST);
    dfm2::opengl::myGlColorDiffuse(dfm2::CColor::Red());
    ::glColor3d(0,0,0);
    dfm2::opengl::DrawCylinder(ps, pe, 0.001);
    dfm2::opengl::DrawAxis(1);
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


