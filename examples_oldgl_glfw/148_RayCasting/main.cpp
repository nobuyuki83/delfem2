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

#include "delfem2/srch_trimesh3_class.h"
#include "delfem2/srch_bruteforce.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/srch_bvh.h"
#include "delfem2/msh_points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mat4.h"
#include "delfem2/thread.h"
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
    const double mMVPd[16],
    const std::vector<delfem2::PointOnSurfaceMesh<double> > &aPointElemSurf,
    const std::vector<double> &aXYZ, // 3d points
    const std::vector<unsigned int> &aTri) {
  const dfm2::CMat4d mMVPd_inv = dfm2::Inverse_Mat4(mMVPd);
  aRGB.resize(nheight * nwidth * 3);
  auto render = [&](int iw, int ih) {
    const std::pair<dfm2::CVec3d, dfm2::CVec3d> ray = dfm2::RayFromInverseMvpMatrix(
        mMVPd_inv.data(), iw, ih, nwidth, nheight);
    const delfem2::PointOnSurfaceMesh<double> &pes = aPointElemSurf[ih * nwidth + iw];
    if (pes.itri == UINT_MAX) {
      aRGB[(ih * nwidth + iw) * 3 + 0] = 200;
      aRGB[(ih * nwidth + iw) * 3 + 1] = 255;
      aRGB[(ih * nwidth + iw) * 3 + 2] = 255;
    } else {
      const unsigned int itri = pes.itri;
      assert(itri < aTri.size() / 3);
      dfm2::CVec3d n = dfm2::Normal_TriInMeshTri3(itri, aXYZ.data(), aTri.data());
      const double dot = n.normalized().dot(ray.second);
      aRGB[(ih * nwidth + iw) * 3 + 0] = static_cast<unsigned char>(-dot * 255);
      aRGB[(ih * nwidth + iw) * 3 + 1] = static_cast<unsigned char>(-dot * 255);
      aRGB[(ih * nwidth + iw) * 3 + 2] = static_cast<unsigned char>(-dot * 255);
    }
  };

  delfem2::parallel_for(nwidth, nheight, render);
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

  std::vector<dfm2::CNodeBVH2> bvh_nodes;
  std::vector<dfm2::CBV3_Sphere<double>> bvh_volumes;
  delfem2::ConstructBVHTriangleMeshMortonCode(
      bvh_nodes, bvh_volumes,
      vec_xyz, vec_tri);

  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    tex.width = 256;
    tex.height = 256;
    tex.channels = 3;
    tex.pixel_color.resize(tex.width * tex.height * tex.channels);
  }

  dfm2::glfw::CViewer3 viewer(2);
  viewer.width = 400;
  viewer.height = 400;
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
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
      const dfm2::CMat4f mP = viewer.GetProjectionMatrix();
      const dfm2::CMat4f mMV = viewer.GetModelViewMatrix();
      const dfm2::CMat4d mMVP = (mP * mMV).cast<double>();
      std::vector<delfem2::PointOnSurfaceMeshd> vec_point_on_tri;
      Intersection_ImageRay_TriMesh3(
          vec_point_on_tri,
          tex.height, tex.width, mMVP.data(),
          bvh_nodes, bvh_volumes, vec_xyz, vec_tri);
      ShadingImageRayLambertian(
          tex.pixel_color,
          tex.height, tex.width, mMVP.data(),
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


