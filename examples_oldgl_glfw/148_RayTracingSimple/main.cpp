/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/tex.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/srchuni_v3.h"
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"
#include "delfem2/mat4.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <algorithm>

namespace dfm2 = delfem2;

// ----------------------------------------

void ShadingImageRayLambertian(
    std::vector<unsigned char>& aRGB,
    unsigned int nheight,
    unsigned int nwidth,
    const float mMVPf[16],
    const std::vector< delfem2::CPtElm2<double> >& aPointElemSurf,
    const std::vector<double>& aXYZ, // 3d points
    const std::vector<unsigned int>& aTri)
{
  double mMVPd[16]; for(int i=0;i<16;++i){ mMVPd[i] = mMVPf[i]; }
  double mMVPd_inv[16]; dfm2::Inverse_Mat4(mMVPd_inv,mMVPd);
  aRGB.resize(nheight*nwidth*3);
  for(unsigned int ih=0;ih<nheight;++ih){
    for(unsigned int iw=0;iw<nwidth;++iw) {
      const double ps[4] = { -1. + (2./nwidth)*(iw+0.5), -1. + (2./nheight)*(ih+0.5), -1., 1. };
      const double pe[4] = { -1. + (2./nwidth)*(iw+0.5), -1. + (2./nheight)*(ih+0.5), +1., 1. };
      double qs[3]; dfm2::Vec3_Vec3Mat4_AffineProjection(qs, ps,mMVPd_inv);
      double qe[3]; dfm2::Vec3_Vec3Mat4_AffineProjection(qe, pe,mMVPd_inv);
      const dfm2::CVec3d src1(qs);
      const dfm2::CVec3d dir1 = dfm2::CVec3d(qe) - src1;
      //
      const delfem2::CPtElm2<double>& pes = aPointElemSurf[ih*nwidth+iw];
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
        dfm2::CVec3d udir1 = dir1.Normalize();
        const double dot = n[0] * udir1[0] + n[1] * udir1[1] + n[2] * udir1[2];
        aRGB[(ih * nwidth + iw) * 3 + 0] = static_cast<unsigned char>(-dot * 255);
        aRGB[(ih * nwidth + iw) * 3 + 1] = static_cast<unsigned char>(-dot * 255);
        aRGB[(ih * nwidth + iw) * 3 + 2] = static_cast<unsigned char>(-dot * 255);
      }
    }
  }
}

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ; // 3d points
  std::vector<unsigned int> aTri;

  { // load input mesh
    delfem2::Read_Ply(std::string(PATH_INPUT_DIR) + "/bunny_2k.ply",
                      aXYZ, aTri);
    dfm2::Normalize_Points3(aXYZ, 2.0);
  }

  std::vector<dfm2::CNodeBVH2> aNodeBVH;
  std::vector<dfm2::CBV3_Sphere<double>> aAABB;
  {
    std::vector<double> aCent;
    dfm2::CentsMaxRad_MeshTri3(aCent,
        aXYZ,aTri);
    double min_xyz[3], max_xyz[3];
    delfem2::BoundingBox3_Points3(min_xyz,max_xyz,
        aCent.data(), aCent.size()/3);
    std::vector<unsigned int> aSortedId;
    std::vector<std::uint32_t> aSortedMc;
    dfm2::SortedMortenCode_Points3(aSortedId,aSortedMc,
                                   aCent,min_xyz,max_xyz);
    dfm2::BVHTopology_Morton(aNodeBVH,
                             aSortedId,aSortedMc);
#ifndef NDEBUG
    dfm2::Check_MortonCode_Sort(aSortedId, aSortedMc, aCent, min_xyz, max_xyz);
    dfm2::Check_MortonCode_RangeSplit(aSortedMc);
#endif
    dfm2::CLeafVolumeMaker_Mesh<dfm2::CBV3d_Sphere,double> lvm(
        1.0e-10,
        aXYZ.data(), aXYZ.size()/3,
        aTri.data(), aTri.size()/3, 3);
    dfm2::BVH_BuildBVHGeometry(aAABB,
        0, aNodeBVH,
        lvm);
#ifndef NDEBUG
      dfm2::Check_BVH(aNodeBVH,aCent.size()/3);
#endif
  }

  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    tex.w = 256;
    tex.h = 256;
    tex.aRGB.resize(tex.w*tex.h*3);
  }

  dfm2::opengl::CViewer_GLFW viewer;
  viewer.width = 400;
  viewer.height = 400;
  viewer.camera.view_height = 2;
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;

  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();

  tex.InitGL();

  while (!glfwWindowShouldClose(viewer.window))
  {
    for(unsigned int i=0;i<10;++i) {
      viewer.DrawBegin_oldGL();
      ::glColor3d(0.8,0.8,0.8);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
      ::glColor3d(0.0,0.0,0.0);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aTri);
      viewer.SwapBuffers();
      glfwPollEvents();
    }
    for(unsigned int i=0;i<10;++i) {
      float mMVP[16];
      {
        float mMV[16], mP[16];
        {
          int width0, height0;
          glfwGetFramebufferSize(viewer.window, &width0, &height0);
          viewer.camera.Mat4_MVP_OpenGL(mMV, mP, float(width0)/float(height0));
        }
        dfm2::MatMat4(mMVP, mMV, mP);
      }
      std::vector< delfem2::CPtElm2d > aPointElemSurf;
      Intersection_ImageRay_TriMesh3(aPointElemSurf,
           tex.h,tex.w, mMVP,
           aNodeBVH,aAABB,aXYZ,aTri);
      ShadingImageRayLambertian(tex.aRGB,
          tex.h, tex.w, mMVP,
          aPointElemSurf, aXYZ, aTri);
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


