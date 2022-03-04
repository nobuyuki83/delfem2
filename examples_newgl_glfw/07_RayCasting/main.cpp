/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#else
#  include <glad/glad.h>
#endif
#if defined(_MSC_VER)
#  include <windows.h>
#endif
#include <GLFW/glfw3.h>

#include "delfem2/srch_trimesh3_class.h"
#include "delfem2/srch_bruteforce.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_normal.h"
#include "delfem2/srch_bvh.h"
#include "delfem2/mat4.h"
#include "delfem2/thread.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/drawer_mshtex.h"
#include "delfem2/opengl/tex.h"

namespace dfm2 = delfem2;

// ----------------------------------------

void ShadingImageRayLambertian(
    std::vector<unsigned char> &vec_rgb,
    unsigned int height,
    unsigned int width,
    const double mat_mvp_rowmajor[16],
    const std::vector<delfem2::PointOnSurfaceMesh<double> > &aPointElemSurf,
    const std::vector<double> &vec_xyz, // 3d points
    const std::vector<unsigned int> &vec_tri,
    bool is_parallel) {
  const dfm2::CMat4d mat_mvp_colmajor_inverse = delfem2::Inverse_Mat4(mat_mvp_rowmajor);
  vec_rgb.resize(height * width * 3);
  auto func0 = [&](int ih, int iw) {
    const std::pair<dfm2::CVec3d, dfm2::CVec3d> ray = delfem2::RayFromInverseMvpMatrix(
        mat_mvp_colmajor_inverse.data(), iw, ih, width, height);
    const delfem2::PointOnSurfaceMesh<double> &pes = aPointElemSurf[ih * width + iw];
    if (pes.itri == UINT_MAX) {
      vec_rgb[(ih * width + iw) * 3 + 0] = 200;
      vec_rgb[(ih * width + iw) * 3 + 1] = 255;
      vec_rgb[(ih * width + iw) * 3 + 2] = 255;
    } else {
      const unsigned int itri = pes.itri;
      assert(itri < vec_tri.size() / 3);
      dfm2::CVec3d n = delfem2::Normal_TriInMeshTri3(itri, vec_xyz.data(), vec_tri.data());
      const double dot = n.normalized().dot(ray.second);
      vec_rgb[(ih * width + iw) * 3 + 0] = static_cast<unsigned char>(-dot * 255);
      vec_rgb[(ih * width + iw) * 3 + 1] = static_cast<unsigned char>(-dot * 255);
      vec_rgb[(ih * width + iw) * 3 + 2] = static_cast<unsigned char>(-dot * 255);
    }
  };
  if (is_parallel) {
    delfem2::parallel_for(width, height, func0);
  } else {
    for (unsigned int ih = 0; ih < height; ++ih) {
      for (unsigned int iw = 0; iw < width; ++iw) {
        func0(ih, iw);
      }
    }
  }
}

struct Data {
  std::vector<double> aXYZ; // 3d points
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CNodeBVH2> aNodeBVH;
  std::vector<dfm2::CBV3_Sphere<double>> aAABB;
  //
  dfm2::opengl::Drawer_RectangleTex drawer;
  dfm2::glfw::CViewer3 viewer;
  dfm2::opengl::CTexRGB_Rect2D tex;
#ifdef EMSCRIPTEN
  const bool is_parallel = false;
#else
  const bool is_parallel = true;
#endif
};

void draw(Data *data) {
  {
    delfem2::CMat4f mP = data->viewer.GetProjectionMatrix();
    delfem2::CMat4f mMV = data->viewer.GetModelViewMatrix();
    delfem2::CMat4d mMVP = (mP * mMV).cast<double>();
    std::vector<delfem2::PointOnSurfaceMesh<double> > aPointElemSurf;
    Intersection_ImageRay_TriMesh3(
        aPointElemSurf,
        data->tex.height,
        data->tex.width,
        mMVP.data(),
        data->aNodeBVH, data->aAABB,
        data->aXYZ, data->aTri,
        data->is_parallel);
    ShadingImageRayLambertian(
        data->tex.pixel_color,
        data->tex.height, data->tex.width, mMVP.data(),
        aPointElemSurf, data->aXYZ, data->aTri,
        data->is_parallel);
  }
  data->tex.InitGL();
  //
  ::glfwMakeContextCurrent(data->viewer.window);
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//  ::glEnable(GL_DEPTH_TEST);
//  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);
  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0); // activate the texture unit first before binding texture
  glBindTexture(GL_TEXTURE_2D, data->tex.id_tex);
  data->drawer.Draw(dfm2::CMat4f::Identity().data(),
                    dfm2::CMat4f::Identity().data());
  data->viewer.SwapBuffers();
  glfwPollEvents();
}

int main() {
  Data data;

  { // load input mesh
    std::string path = std::filesystem::path(PATH_INPUT_DIR) / "bunny_2k.ply";
    delfem2::Read_Ply(
        data.aXYZ, data.aTri,
        path);
    std::cout << data.aTri.size() / 3 << std::endl;
    dfm2::Normalize_Points3(data.aXYZ, 2.0);
  }
  dfm2::ConstructBVHTriangleMeshMortonCode(
      data.aNodeBVH, data.aAABB,
      data.aXYZ, data.aTri);

  {
    data.tex.width = 256;
    data.tex.height = 256;
    data.tex.channels = 3;
    data.tex.pixel_color.resize(data.tex.width * data.tex.height * data.tex.channels);
    data.viewer.projection = std::make_unique<delfem2::Projection_LookOriginFromZplus>(2, false);
    data.viewer.width = 400;
    data.viewer.height = 400;
  }
  // -----------------------
  // opengl starts from here
  delfem2::glfw::InitGLNew();
  data.viewer.OpenWindow();
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  data.tex.InitGL();
  data.drawer.InitGL();

#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, &data, 0, 1);
#else
  while (!glfwWindowShouldClose(data.viewer.window)) { draw(&data); }
#endif

  glfwDestroyWindow(data.viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


