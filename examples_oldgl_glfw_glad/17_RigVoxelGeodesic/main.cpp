/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should be before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "delfem2/geo_plane.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/gridvoxel.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/r2tglo.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

void Draw_CGrid3(
    const dfm2::CGrid3<int> &grid,
    [[maybe_unused]] const std::vector<std::pair<double, dfm2::CColor> > &colorMap,
    [[maybe_unused]] double thresh) {
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
          dfm2::opengl::myGlMaterialDiffuse(dfm2::CColor::Red());
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

int main() {
  std::vector<double> vtx_xyz_ini;
  std::vector<unsigned int> tri_vtx;

  dfm2::Read_Ply(
      vtx_xyz_ini, tri_vtx,
      std::filesystem::path(PATH_INPUT_DIR) / "arm_16k.ply");
  dfm2::Normalize_Points3(vtx_xyz_ini, 1.0);

  std::vector<dfm2::CRigBone> bones;
  { // 0
    dfm2::CRigBone b;
    b.ibone_parent = -1;
    b.transRelative[0] = -0.0;
    b.transRelative[1] = -0.0;
    b.transRelative[2] = -0.0;
    bones.push_back(b);
  }
  { // 1
    dfm2::CRigBone b;
    b.ibone_parent = 0;
    b.transRelative[0] = +0.01;
    b.transRelative[1] = +0.20;
    b.transRelative[2] = +0.01;
    bones.push_back(b);
  }
  { // 2
    dfm2::CRigBone b;
    b.ibone_parent = 1;
    b.transRelative[0] = +0.20;
    b.transRelative[1] = +0.01;
    b.transRelative[2] = +0.01;
    bones.push_back(b);
  }
  { // 3
    dfm2::CRigBone b;
    b.ibone_parent = 2;
    b.transRelative[0] = +0.14;
    b.transRelative[1] = +0.16;
    b.transRelative[2] = -0.21;
    bones.push_back(b);
  }
  { // 4
    dfm2::CRigBone b;
    b.ibone_parent = 1;
    b.transRelative[0] = -0.20;
    b.transRelative[1] = +0.01;
    b.transRelative[2] = +0.01;
    bones.push_back(b);
  }
  { // 5
    dfm2::CRigBone b;
    b.ibone_parent = 4;
    b.transRelative[0] = -0.14;
    b.transRelative[1] = +0.10;
    b.transRelative[2] = -0.21;
    bones.push_back(b);
  }
  { // 6
    dfm2::CRigBone b;
    b.ibone_parent = 0;
    b.transRelative[0] = +0.15;
    b.transRelative[1] = -0.2;
    b.transRelative[2] = +0.01;
    bones.push_back(b);
  }
  { // 7
    dfm2::CRigBone b;
    b.ibone_parent = 6;
    b.transRelative[0] = +0.01;
    b.transRelative[1] = -0.2;
    b.transRelative[2] = +0.15;
    bones.push_back(b);
  }
  { // 8
    dfm2::CRigBone b;
    b.ibone_parent = 0;
    b.transRelative[0] = -0.15;
    b.transRelative[1] = -0.2;
    b.transRelative[2] = +0.01;
    bones.push_back(b);
  }
  { // 9
    dfm2::CRigBone b;
    b.ibone_parent = 8;
    b.transRelative[0] = -0.01;
    b.transRelative[1] = -0.2;
    b.transRelative[2] = +0.15;
    bones.push_back(b);
  }
  dfm2::UpdateBoneRotTrans(bones);
  for (auto & ib : bones) {
    dfm2::Mat4_AffineTranslation(
        ib.invBindMat,
        -ib.affmat3Global[0 * 4 + 3],
        -ib.affmat3Global[1 * 4 + 3],
        -ib.affmat3Global[2 * 4 + 3]);
  }

  // ---------------------------------------

  dfm2::opengl::CRender2Tex_DrawOldGL_BOX sampler_box;
  {
    unsigned int ndiv = 81;
    sampler_box.Initialize(ndiv, ndiv, ndiv, 1.0 / ndiv);
  }

  // ---------------------------------------
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
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
    dfm2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz_ini, tri_vtx);
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

  std::vector<double> aW;
  {
    const unsigned int np = vtx_xyz_ini.size() / 3;
    const unsigned int nb = bones.size();
    aW.resize(np * nb, 0.0);
    for (unsigned ibone = 0; ibone < nb; ++ibone) { // distance from line
      const int ibp = bones[ibone].ibone_parent;
      if (ibp == -1) { continue; }
      const dfm2::CVec3d ps(-bones[ibp].invBindMat[3],
                            -bones[ibp].invBindMat[7],
                            -bones[ibp].invBindMat[11]);
      const dfm2::CVec3d pe(-bones[ibone].invBindMat[3],
                            -bones[ibone].invBindMat[7],
                            -bones[ibone].invBindMat[11]);
//      std::cout << ps << " " << pe << std::endl;
      std::vector<std::pair<unsigned int, double> > aIdvoxDist;
      {
        std::vector<unsigned int> aIndvox;
        Intersection_VoxelGrid_LinSeg(aIndvox,
                                      grid, ps, pe);
//        std::cout << "aIndvox.size(): " << aIndvox.size() << std::endl;
        const unsigned int nx = grid.ndivx;
        const unsigned int ny = grid.ndivy;
        for (unsigned int ivox0: aIndvox) {
          if (grid.aVal[ivox0] == 0) { continue; }
          const unsigned int iz0 = ivox0 / (ny * nx);
          const unsigned int iy0 = (ivox0 - iz0 * ny * nx) / nx;
          const unsigned int ix0 = ivox0 - iz0 * ny * nx - iy0 * nx;
          dfm2::CVec3d p0(ix0 + 0.5, iy0 + 0.5, iz0 + 0.5);
          dfm2::CVec3d p1;
          dfm2::Vec3_Mat4Vec3_Affine(p1.p, grid.am.mat, p0.p);
          dfm2::CVec3d p2 = dfm2::Nearest_Edge_Point(p1, ps, pe);
          double dist0 = (p2 - p1).norm();
          aIdvoxDist.emplace_back(ivox0, dist0);
        }
      }
      std::vector<double> aDist;
      VoxelGeodesic(aDist,
                    aIdvoxDist, sampler_box.edgeLen(), grid);
      const unsigned int nx = grid.ndivx;
      const unsigned int ny = grid.ndivy;
      const unsigned int nz = grid.ndivz;
      const dfm2::CMat4d &ami = grid.am.Inverse();
      for (unsigned int ip = 0; ip < np; ip++) {
        dfm2::CVec3d p0;
        dfm2::Vec3_Mat4Vec3_Affine(p0.p, ami.mat, vtx_xyz_ini.data() + ip * 3);
        int ix0 = (int) floor(p0.x - 0.5);
        int iy0 = (int) floor(p0.y - 0.5);
        int iz0 = (int) floor(p0.z - 0.5);
        double dist0 = 0;
        double sumw = 0.0;
        for (unsigned int i = 0; i < 8; ++i) {
          int ix = ix0;
          if (i % 2 == 1) { ix += 1; }
          int iy = iy0;
          if ((i / 2) % 2 == 1) { iy += 1; }
          int iz = iz0;
          if ((i / 4) % 2 == 1) { iz += 1; }
          assert(ix >= -1 && ix <= (int)nx);
          assert(iy >= -1 && iy <= (int)ny);
          assert(iz >= -1 && iz <= (int)nz);
          unsigned int ivox = iz * ny * nx + iy * nx + ix;
          if (grid.aVal[ivox] == 0) continue;
          double len = (p0.x - ix) * (p0.x - ix) + (p0.y - iy) * (p0.y - iy) + (p0.z - iz) * (p0.z - iz);
          len = sqrt(len);
          double w0 = 1.0 / (len + 1.0e-3);
//a          std::cout << ibone << " " << ip << " " << aDist[ivox] << std::endl;
          dist0 += aDist[ivox];
          sumw += w0;
        }
        assert(sumw != 0.0);
        dist0 = dist0 / sumw;
        double tmp = dist0 * dist0 + 1.0e-10;
        aW[ip * nb + ibp] = 1.0 / (tmp * tmp);
      }
    }
  } // nbone
  {
    unsigned int nb = bones.size();
    for (unsigned int ip = 0; ip < vtx_xyz_ini.size() / 3; ++ip) {
      double w0 = 0.0;
      for (unsigned int ib = 0; ib < nb; ++ib) {
        w0 += aW[ip * nb + ib];
      }
      w0 = 1.0 / w0;
      for (unsigned int ib = 0; ib < nb; ++ib) {
        aW[ip * nb + ib] *= w0;
//        std::cout << "   " << aW[ip*nb+ib] << std::endl;
      }
    }
  }

  std::vector<double> aXYZ1 = vtx_xyz_ini;

  std::random_device rd;
  std::mt19937 reng(rd());
  std::uniform_int_distribution<unsigned int> dist0(0, grid.aVal.size() - 1);
  std::uniform_real_distribution<double> dist1(-0.5 * sampler_box.nDivX() * sampler_box.edgeLen(),
                                               +0.5 * sampler_box.nDivX() * sampler_box.edgeLen());

  std::vector<std::pair<double, dfm2::CColor> > colorMap;
  { // set color map
    double el = sampler_box.edgeLen();
    const unsigned int nx = sampler_box.nDivX();
    const unsigned int ny = sampler_box.nDivY();
    const unsigned int nz = sampler_box.nDivZ();
    const double len = el * sqrt(nx * nx + ny * ny + nz * nz);
    dfm2::ColorMap_BlueCyanGreenYellowRed(colorMap, 0, len);
  }


  // -------
  unsigned int iframe = 0;
  while (true) {
    iframe = 0;
    for (; iframe < 100; ++iframe) { // draw the result
      bones[1].SetRotationBryant(0, sin(iframe * 0.01), 0.0);
      dfm2::UpdateBoneRotTrans(bones);
      dfm2::Skinning_LBS(aXYZ1, vtx_xyz_ini, bones, aW);
      // --------
      viewer.DrawBegin_oldGL();
//      Draw_CGrid3(grid,colorMap,iframe*0.1);
      dfm2::opengl::myGlColorDiffuse(dfm2::CColor::Gray(0.9));
      ::glEnable(GL_LIGHTING);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), tri_vtx.data(), tri_vtx.size() / 3);
      ::glDisable(GL_DEPTH_TEST);
      ::glDisable(GL_LIGHTING);
      dfm2::opengl::DrawBoneCurrent_Line(bones, -1, 0.02);
      dfm2::opengl::myGlColorDiffuse(dfm2::CColor::Red());
      ::glColor3d(0, 0, 0);
      viewer.SwapBuffers();
      glfwPollEvents();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


