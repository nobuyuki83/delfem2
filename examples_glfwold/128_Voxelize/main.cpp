/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/gridvoxel.h"
//
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/render2tex_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------


void Draw_CGrid3
(const dfm2::CGrid3<int>& grid)
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
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::myGlMaterialDiffuse(dfm2::CColor::Gray(0.8));
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
  // --------
  // end transformation
  ::glPopMatrix();
}


// ------------------------------------------------------

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;

  dfm2::Read_Obj(
      std::string(PATH_INPUT_DIR)+"/bunny_1k.obj",
      aXYZ,aTri);
  dfm2::Normalize_Points3(aXYZ,4.0);
  // ---------------------------------------

  dfm2::opengl::CRender2Tex_DrawOldGL_BOX sampler_box;
  sampler_box.Initialize(64, 64, 64, 0.075);

  // ---------------------------------------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 3.0;
  viewer.nav.camera.camera_rot_mode = dfm2::CAMERA_ROT_TBALL;
//  viewer.nav.camera.Rot_Camera(+0.2, -0.2);
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);
  // ------------
  for(auto& smplr: sampler_box.aSampler){
    smplr.InitGL(); // move the sampled image to a texture
    smplr.Start();
    ::glClearColor(1.0, 1.0, 1.0, 1.0 );
    ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
    smplr.End();
    smplr.GetDepth();
  }
  
  dfm2::CGrid3<int> grid;
  {
    const unsigned int nx = sampler_box.nDivX();
    const unsigned int ny = sampler_box.nDivY();
    const unsigned int nz = sampler_box.nDivZ();
    const double el = sampler_box.edgeLen();
    grid.am.SetScale(el, el, el);
    grid.am.mat[ 3] = -(nx*el*0.5);
    grid.am.mat[ 7] = -(ny*el*0.5);
    grid.am.mat[11] = -(nz*el*0.5);
    grid.Initialize(nx,ny,nz, 1);
    for(unsigned int iy=0;iy<ny;++iy){
      for(unsigned int iz=0;iz<nz;++iz){
        double d0 = sampler_box.aSampler[0].aZ[iz*ny+iy];
        const unsigned int nd = floor(-d0/el+1.0e-5);
        for(unsigned int id=0;id<nd;id++){
          const unsigned int ix0 = nx-1-id;
          const unsigned int iy0 = iy;
          const unsigned int iz0 = iz;
          grid.aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
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
          grid.aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
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
          grid.aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
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
          grid.aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
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
          grid.aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
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
          grid.aVal[iz0*ny*nx+iy0*nx+ix0] = 0;
        }
      }
    }
  }

  while (true)
  {
    viewer.DrawBegin_oldGL();
    sampler_box.Draw();
    Draw_CGrid3(grid);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


