/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should be before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/r2tglo.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/gridvoxel.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

void Draw_CGrid3
    (const dfm2::CGrid3<int> &grid) {
  { // set-up transformation
    const dfm2::CMat4d &am = grid.am;
    dfm2::CMat4d amt = am.transpose();
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
    for (unsigned int iz = 0; iz < nz; ++iz) {
      for (unsigned int iy = 0; iy < ny; ++iy) {
        for (unsigned int ix = 0; ix < nx; ++ix) {
          if (grid.aVal[iz * ny * nx + iy * nx + ix] == 0) { continue; }
          const double pmin[3] = {(double) ix, (double) iy, (double) iz};
          const double pmax[3] = {(double) ix + 1., (double) iy + 1., (double) iz + 1.};
          dfm2::opengl::DrawBox3_Face(pmin, pmax);
        }
      }
    }
  }
  // --------
  // end transformation
  ::glPopMatrix();
}


// ------------------------------------------------------

int main(int argc, char *argv[]) {
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;

  dfm2::Read_Obj(
      std::string(PATH_INPUT_DIR) + "/bunny_1k.obj",
      aXYZ, aTri);
  dfm2::Normalize_Points3(aXYZ, 4.0);
  // ---------------------------------------

  dfm2::opengl::CRender2Tex_DrawOldGL_BOX sampler_box;
  sampler_box.Initialize(128, 128, 128, 0.04);

  // ---------------------------------------
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 3.0;
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
//  viewer.camera.Rot_Camera(+0.2, -0.2);
  if (!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);
  // ------------
  for (auto &smplr: sampler_box.aSampler) {
    smplr.InitGL(); // move the sampled image to a texture
    smplr.Start();
    dfm2::opengl::SetView(smplr);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri);
    smplr.End();
  }

  dfm2::CGrid3<int> grid;
  {
    const unsigned int nx = sampler_box.nDivX();
    const unsigned int ny = sampler_box.nDivY();
    const unsigned int nz = sampler_box.nDivZ();
    const double el = sampler_box.edgeLen();
    grid.am.SetScale(el, el, el);
    grid.am.mat[3] = -(nx * el * 0.5);
    grid.am.mat[7] = -(ny * el * 0.5);
    grid.am.mat[11] = -(nz * el * 0.5);
    grid.Initialize(nx, ny, nz, 1);
  }
  CarveVoxelByDepth(grid.aVal,
                    sampler_box);

  while (true) {
    for (int itr = 0; itr < 5; ++itr) { // dilation
      dfm2::Grid3Voxel_Dilation(grid);
      viewer.DrawBegin_oldGL();
      sampler_box.Draw();
      Draw_CGrid3(grid);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if (glfwWindowShouldClose(viewer.window)) goto EXIT;
    }
    for (int itr = 0; itr < 5; ++itr) { // erosion
      dfm2::Grid3Voxel_Erosion(grid);
      viewer.DrawBegin_oldGL();
      sampler_box.Draw();
      Draw_CGrid3(grid);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if (glfwWindowShouldClose(viewer.window)) goto EXIT;
    }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


