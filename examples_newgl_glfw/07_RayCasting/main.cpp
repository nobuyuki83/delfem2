/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#if defined(_MSC_VER)
#  include <windows.h>
#endif
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/srchuni_v3.h"
#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshmisc.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"
#include "delfem2/mat4.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/drawer_mshtex.h"
#include "delfem2/opengl/tex.h"

namespace dfm2 = delfem2;

// ----------------------------------------

void ShadingImageRayLambertian(
    std::vector<unsigned char>& vec_rgb,
    unsigned int height,
    unsigned int width,
    const float mat_mvp_rowmajor_float[16],
    const std::vector< delfem2::PointOnSurfaceMesh<double> >& aPointElemSurf,
    const std::vector<double>& vec_xyz, // 3d points
    const std::vector<unsigned int>& vec_tri)
{
  double mat_mvp_rowmajor_double[16];
  for(int i=0;i<16;++i){ mat_mvp_rowmajor_double[i] = mat_mvp_rowmajor_float[i]; }
  double mat_mvp_colmajor_inverse[16];
  dfm2::Inverse_Mat4(mat_mvp_colmajor_inverse,mat_mvp_rowmajor_double);
  vec_rgb.resize(height*width*3);
  for(unsigned int ih=0;ih<height;++ih){
    for(unsigned int iw=0;iw<width;++iw) {
      const double ps[4] = { -1. + (2./width)*(iw+0.5), -1. + (2./height)*(ih+0.5), +1., 1. };
      const double pe[4] = { -1. + (2./width)*(iw+0.5), -1. + (2./height)*(ih+0.5), -1., 1. };
      double qs[3]; dfm2::Vec3_Mat4Vec3_Homography(qs, mat_mvp_colmajor_inverse, ps);
      double qe[3]; dfm2::Vec3_Mat4Vec3_Homography(qe, mat_mvp_colmajor_inverse, pe);
      const dfm2::CVec3d src1(qs);
      const dfm2::CVec3d dir1 = dfm2::CVec3d(qe) - src1;
      //
      const delfem2::PointOnSurfaceMesh<double>& pes = aPointElemSurf[ih*width+iw];
      if (pes.itri == UINT_MAX) {
        vec_rgb[(ih * width + iw) * 3 + 0] = 200;
        vec_rgb[(ih * width + iw) * 3 + 1] = 255;
        vec_rgb[(ih * width + iw) * 3 + 2] = 255;
      } else {
        const unsigned int itri = pes.itri;
        assert(itri < vec_tri.size() / 3);
        double n[3], area;
        delfem2::UnitNormalAreaTri3(
            n, area,
            vec_xyz.data() + vec_tri[itri * 3 + 0] * 3,
            vec_xyz.data() + vec_tri[itri * 3 + 1] * 3,
            vec_xyz.data() + vec_tri[itri * 3 + 2] * 3);
        dfm2::CVec3d udir1 = dir1.normalized();
        const double dot = n[0] * udir1[0] + n[1] * udir1[1] + n[2] * udir1[2];
        vec_rgb[(ih * width + iw) * 3 + 0] = static_cast<unsigned char>(-dot * 255);
        vec_rgb[(ih * width + iw) * 3 + 1] = static_cast<unsigned char>(-dot * 255);
        vec_rgb[(ih * width + iw) * 3 + 2] = static_cast<unsigned char>(-dot * 255);
      }
    }
  }
}

int main()
{
  std::vector<double> aXYZ; // 3d points
  std::vector<unsigned int> aTri;

  { // load input mesh
    delfem2::Read_Ply(
        aXYZ, aTri,
        std::filesystem::path(PATH_INPUT_DIR) / "bunny_2k.ply");
    dfm2::Normalize_Points3(aXYZ, 2.0);
  }

  std::vector<dfm2::CNodeBVH2> aNodeBVH;
  std::vector<dfm2::CBV3_Sphere<double>> aAABB;
  dfm2::ConstructBVHTriangleMeshMortonCode(
      aNodeBVH, aAABB,
      aXYZ, aTri);
  
  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    tex.width = 256;
    tex.height = 256;
    tex.channels = 3;
    tex.pixel_color.resize(tex.width*tex.height*tex.channels);
  }

  dfm2::opengl::Drawer_RectangleTex drawer;

  dfm2::glfw::CViewer3 viewer(2);
  viewer.width = 400;
  viewer.height = 400;

  // -----------------------
  // opengl starts from here
  delfem2::glfw::InitGLNew();
  viewer.OpenWindow();
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
  tex.InitGL();
  drawer.InitGL();

  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      delfem2::CMat4f mP = viewer.GetProjectionMatrix();
      delfem2::CMat4f mMV = viewer.GetModelViewMatrix();
      delfem2::CMat4f mMVP = mP * mMV;
      std::vector< delfem2::PointOnSurfaceMesh<double> > aPointElemSurf;
      Intersection_ImageRay_TriMesh3(aPointElemSurf,
           tex.height,tex.width, mMVP.data(),
           aNodeBVH,aAABB,aXYZ,aTri);
      ShadingImageRayLambertian(tex.pixel_color,
          tex.height, tex.width, mMVP.data(),
          aPointElemSurf, aXYZ, aTri);
    }
    tex.InitGL();
    //
    ::glfwMakeContextCurrent(viewer.window);
    ::glClearColor(0.8, 1.0, 1.0, 1.0);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//  ::glEnable(GL_DEPTH_TEST);
//  ::glDepthFunc(GL_LESS);
    ::glEnable(GL_POLYGON_OFFSET_FILL );
    ::glPolygonOffset( 1.1f, 4.0f );

    glEnable(GL_TEXTURE_2D);
    glActiveTexture(GL_TEXTURE0); // activate the texture unit first before binding texture
    glBindTexture(GL_TEXTURE_2D , tex.id_tex);
    drawer.Draw(dfm2::CMat4f::Identity().data(),
                dfm2::CMat4f::Identity().data());
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


