/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/srchuni_v3.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"
#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mat4.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/tex.h"

namespace dfm2 = delfem2;

// ----------------------------------------

void ShadingImageRayLambertian(
    std::vector<unsigned char> &aRGB,
    unsigned int nheight,
    unsigned int nwidth,
    const float mMVPf[16],
    const std::vector<delfem2::PointOnSurfaceMesh<double> > &aPointElemSurf,
    const std::vector<double> &aXYZ, // 3d points
    const std::vector<unsigned int> &aTri) {
  double mMVPd[16];
  for (int i = 0; i < 16; ++i) { mMVPd[i] = mMVPf[i]; }
  double mMVPd_inv[16];
  dfm2::Inverse_Mat4(mMVPd_inv, mMVPd);
  aRGB.resize(nheight * nwidth * 3);
  for (unsigned int ih = 0; ih < nheight; ++ih) {
    for (unsigned int iw = 0; iw < nwidth; ++iw) {
      const double ps[4] = {-1. + (2. / nwidth) * (iw + 0.5), -1. + (2. / nheight) * (ih + 0.5), -1., 1.};
      const double pe[4] = {-1. + (2. / nwidth) * (iw + 0.5), -1. + (2. / nheight) * (ih + 0.5), +1., 1.};
      double qs[3];
      dfm2::Vec3_Vec3Mat4_AffineProjection(qs, ps, mMVPd_inv);
      double qe[3];
      dfm2::Vec3_Vec3Mat4_AffineProjection(qe, pe, mMVPd_inv);
      const dfm2::CVec3d src1(qs);
      const dfm2::CVec3d dir1 = dfm2::CVec3d(qe) - src1;
      //
      const delfem2::PointOnSurfaceMesh<double> &pes = aPointElemSurf[ih * nwidth + iw];
      if (pes.itri == UINT_MAX) {
        aRGB[(ih * nwidth + iw) * 3 + 0] = 200;
        aRGB[(ih * nwidth + iw) * 3 + 1] = 255;
        aRGB[(ih * nwidth + iw) * 3 + 2] = 255;
      } else {
        const unsigned int itri = pes.itri;
        assert(itri < aTri.size() / 3);
        double n[3], area;
        delfem2::UnitNormalAreaTri3(
            n, area,
            aXYZ.data() + aTri[itri * 3 + 0] * 3,
            aXYZ.data() + aTri[itri * 3 + 1] * 3,
            aXYZ.data() + aTri[itri * 3 + 2] * 3);
        dfm2::CVec3d udir1 = dir1.normalized();
        const double dot = n[0] * udir1[0] + n[1] * udir1[1] + n[2] * udir1[2];
        aRGB[(ih * nwidth + iw) * 3 + 0] = static_cast<unsigned char>(-dot * 255);
        aRGB[(ih * nwidth + iw) * 3 + 1] = static_cast<unsigned char>(-dot * 255);
        aRGB[(ih * nwidth + iw) * 3 + 2] = static_cast<unsigned char>(-dot * 255);
      }
    }
  }
}

int main() {
  std::vector<double> vec_xyz; // 3d points
  std::vector<unsigned int> vec_tri;

  { // load input mesh
    delfem2::Read_Ply(
        vec_xyz, vec_tri,
        std::filesystem::path(PATH_INPUT_DIR) / "bunny_2k.ply");
    dfm2::Normalize_Points3(vec_xyz, 2.0);
  }

  std::vector<dfm2::CNodeBVH2> vec_node_bvh;
  std::vector<dfm2::CBV3_Sphere<double>> vec_bv;
  delfem2::ConstructBVHTriangleMeshMortonCode(
      vec_node_bvh, vec_bv,
      vec_xyz, vec_tri);

  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    tex.width = 256;
    tex.height = 256;
    tex.channels = 3;
    tex.pixel_color.resize(tex.width * tex.height * tex.channels);
  }

  dfm2::glfw::CViewer3 viewer;
  viewer.width = 400;
  viewer.height = 400;
  viewer.projection.view_height = 2;

  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  tex.InitGL();

  while (!glfwWindowShouldClose(viewer.window)) {
    for (unsigned int i = 0; i < 10; ++i) {
      viewer.DrawBegin_oldGL();
      ::glColor3d(0.8, 0.8, 0.8);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(vec_xyz, vec_tri);
      ::glColor3d(0.0, 0.0, 0.0);
      dfm2::opengl::DrawMeshTri3D_Edge(vec_xyz, vec_tri);
      viewer.SwapBuffers();
      glfwPollEvents();
    }
    for (unsigned int i = 0; i < 10; ++i) {
      float mMVP[16];
      {
        float mMV[16], mP[16];
        {
          int width0, height0;
          glfwGetFramebufferSize(viewer.window, &width0, &height0);
          viewer.projection.Mat4ColumnMajor(mP, float(width0) / float(height0));
          viewer.modelview.Mat4ColumnMajor(mMV);
//          viewer.camera.Mat4_MVP_OpenGL(mMV, mP, );
        }
        dfm2::MatMat4(mMVP, mMV, mP);
      }
      std::vector<delfem2::PointOnSurfaceMeshd> vec_point_on_tri;
      Intersection_ImageRay_TriMesh3(
          vec_point_on_tri,
          tex.height, tex.width, mMVP,
          vec_node_bvh, vec_bv, vec_xyz, vec_tri);
      ShadingImageRayLambertian(
          tex.pixel_color,
          tex.height, tex.width, mMVP,
          vec_point_on_tri, vec_xyz, vec_tri);
      tex.InitGL();
      //
      viewer.DrawBegin_oldGL();
      ::glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
      ::glClear(GL_COLOR_BUFFER_BIT);
      ::glDisable(GL_LIGHTING);
      ::glMatrixMode(GL_PROJECTION);
      ::glLoadIdentity();
      ::glMatrixMode(GL_MODELVIEW);
      ::glLoadIdentity();
      tex.Draw_oldGL();
      viewer.SwapBuffers();
      glfwPollEvents();
    }
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


