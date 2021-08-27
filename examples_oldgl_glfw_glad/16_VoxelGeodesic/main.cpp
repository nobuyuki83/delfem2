/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>
#include <climits>
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
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/r2tglo.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/geoproximity3_v3.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/gridvoxel.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

void Draw_CGrid3(
    const dfm2::CGrid3<int> &grid,
    const std::vector<double> &aDist,
    const std::vector<std::pair<double, dfm2::CColor> > &colorMap,
    double thresh) {
  { // set-up transformation
    const dfm2::CMat4d &am = grid.am;
    dfm2::CMat4d amt = am.transpose();
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    ::glMultMatrixd(amt.mat);
  }
  {
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_NORMALIZE);
    const unsigned int nx = grid.ndivx;
    const unsigned int ny = grid.ndivy;
    const unsigned int nz = grid.ndivz;
    for (unsigned int iz = 0; iz < nz; ++iz) {
      for (unsigned int iy = 0; iy < ny; ++iy) {
        for (unsigned int ix = 0; ix < nx; ++ix) {
          if (grid.aVal[iz * ny * nx + iy * nx + ix] == 0) { continue; }
          const double val = aDist[iz * ny * nx + iy * nx + ix];
          if (val > thresh) { continue; }
          dfm2::CColor c = dfm2::getColor(val, colorMap);
          dfm2::opengl::myGlMaterialDiffuse(c);
          const double pmin[3] = {(double) ix, (double) iy, (double) iz};
          const double pmax[3] = {(double) ix + 1., (double) iy + 1., (double) iz + 1.};
          dfm2::opengl::DrawBox3_Face(pmin, pmax);
        }
      }
    }
  }
  ::glPopMatrix();
}


// ------------------------------------------------------

int main(int argc, char *argv[]) {
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;

  dfm2::Read_Obj(
      std::string(PATH_INPUT_DIR) + "/bunny_1k.obj",
      aXYZ, aTri);
  dfm2::Normalize_Points3(aXYZ, 1.0);
  // ---------------------------------------

  dfm2::opengl::CRender2Tex_DrawOldGL_BOX sampler_box;
  {
    unsigned int ndiv = 80;
    sampler_box.Initialize(ndiv, ndiv, ndiv, 1.0 / ndiv);
  }

  // ---------------------------------------
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.0;
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
    ::glClearColor(1.0, 1.0, 1.0, 1.0);
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
  dfm2::Grid3Voxel_Dilation(grid);

  std::random_device rd;
  std::mt19937 reng(rd());
  std::uniform_int_distribution<unsigned int> dist0(0, grid.aVal.size() - 1);
  std::uniform_real_distribution<double> dist1(-0.5 * sampler_box.nDivX() * sampler_box.edgeLen(),
                                               +0.5 * sampler_box.nDivX() * sampler_box.edgeLen());

  std::vector<std::pair<double, dfm2::CColor> > colorMap;
  {
    double el = sampler_box.edgeLen();
    const unsigned int nx = sampler_box.nDivX();
    const unsigned int ny = sampler_box.nDivY();
    const unsigned int nz = sampler_box.nDivZ();
    const double len = el * sqrt(nx * nx + ny * ny + nz * nz);
    dfm2::ColorMap_BlueCyanGreenYellowRed(colorMap, 0, len);
  }
  // -------
  std::vector<double> aDist;
  while (true) {
    { // distance from a point
      unsigned int ivox0 = UINT_MAX;
      for (int itr = 0; itr < 100; ++itr) {
        ivox0 = dist0(reng);
        if (grid.aVal[ivox0] == 1) { break; }
      }
      if (ivox0 == UINT_MAX) { goto EXIT; }
      std::vector<std::pair<unsigned int, double> > aIdvoxDist;
      aIdvoxDist.emplace_back(ivox0, 0.0);
      VoxelGeodesic(aDist,
                    aIdvoxDist, sampler_box.edgeLen(), grid);
      // ------
      for (int iframe = 0; iframe < 16; ++iframe) { // draw the result
        viewer.DrawBegin_oldGL();
        Draw_CGrid3(grid, aDist, colorMap, iframe * 0.1);
        { // set-up transformation
          const dfm2::CMat4d &am = grid.am;
          dfm2::CMat4d amt = am.transpose();
          ::glMatrixMode(GL_MODELVIEW);
          ::glPushMatrix();
          ::glMultMatrixd(amt.mat);
          ::glDisable(GL_LIGHTING);
          ::glColor3d(1, 0, 0);
          ::glDisable(GL_DEPTH_TEST);
          const unsigned int nx = grid.ndivx;
          const unsigned int ny = grid.ndivy;
          const int iz0 = ivox0 / (ny * nx);
          const int iy0 = (ivox0 - iz0 * ny * nx) / nx;
          const int ix0 = ivox0 - iz0 * ny * nx - iy0 * nx;
          dfm2::opengl::DrawSphereAt(32, 32, 1.0, ix0 + 0.5, iy0 + 0.5, iz0 + 0.5);
          ::glPopMatrix();
        }
        viewer.SwapBuffers();
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
      }
    }
    { // distance from line
      std::vector<std::pair<unsigned int, double> > aIdvoxDist;
      dfm2::CVec3d ps(dist1(reng), dist1(reng), dist1(reng));
      dfm2::CVec3d pe(dist1(reng), dist1(reng), dist1(reng));
      {
        std::vector<unsigned int> aIndvox;
        Intersection_VoxelGrid_LinSeg(aIndvox,
                                      grid, ps, pe);
        const unsigned int nx = grid.ndivx;
        const unsigned int ny = grid.ndivy;
        for (unsigned int ivox0: aIndvox) {
          if (grid.aVal[ivox0] == 0) { continue; }
          const int iz0 = ivox0 / (ny * nx);
          const int iy0 = (ivox0 - iz0 * ny * nx) / nx;
          const int ix0 = ivox0 - iz0 * ny * nx - iy0 * nx;
          dfm2::CVec3d p0(ix0 + 0.5, iy0 + 0.5, iz0 + 0.5);
          dfm2::CVec3d p1;
          dfm2::Vec3_Mat4Vec3_Affine(p1.p, grid.am.mat, p0.p);
          dfm2::CVec3d p2 = dfm2::nearest_LineSeg_Point(p1, ps, pe);
          double dist0 = (p2 - p1).norm();
          aIdvoxDist.emplace_back(ivox0, dist0);
        }
      }
      VoxelGeodesic(aDist,
                    aIdvoxDist, sampler_box.edgeLen(), grid);
      // ------
      for (int iframe = 0; iframe < 16; ++iframe) { // draw the result
        viewer.DrawBegin_oldGL();
        Draw_CGrid3(grid, aDist, colorMap, iframe * 0.1);
        ::glDisable(GL_DEPTH_TEST);
        ::glDisable(GL_LIGHTING);
        dfm2::opengl::myGlColorDiffuse(dfm2::CColor::Red());
        ::glColor3d(0, 0, 0);
        dfm2::opengl::DrawCylinder(ps, pe, 0.005);
        viewer.SwapBuffers();
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
      }
    }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


