/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#include "delfem2/srchuni_v3.h"
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/mat4.h"
// ---------------------------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/tex_gl.h"

namespace dfm2 = delfem2;

// ----------------------------------------


void RayTracing(
    std::vector<unsigned char>& aRGB,
    unsigned int nheight,
    unsigned int nwidth,
    const float mMV[16],
    const std::vector<dfm2::CNodeBVH2>& aNodeBVH,
    const std::vector<dfm2::CBV3_Sphere<double>>& aAABB,
    const std::vector<double>& aXYZ, // 3d points
    const std::vector<unsigned int>& aTri )
{
  double mMVd[16]; for(int i=0;i<16;++i){ mMVd[i] = mMV[i]; }
  double mMVinv[16]; dfm2::Inverse_Mat4(mMVinv,mMVd);
  //
  const double dir0[4] = {0,0,-1,0};
  double dir1[4];
  {
    dfm2::VecMat4(dir1,dir0,mMVinv);
    double l1 = dfm2::Length3(dir1);
    dir1[0] /= l1;
    dir1[1] /= l1;
    dir1[2] /= l1;
  }
  // -----------
  std::vector<unsigned int> aIndElem;
  for(unsigned int ih=0;ih<nheight;++ih){
    for(unsigned int iw=0;iw<nwidth;++iw){
      const double src0[4] = {
          (float)(-2.f + 4.f/(float)nwidth*(iw+0.5)),
          (float)(-2.f + 4.f/(float)nheight*(ih+0.5)),
          (float)(2),
          (float)(1),
      };
      double src1[4]; dfm2::VecMat4(src1,src0,mMVinv);
      aIndElem.resize(0);
      dfm2::BVH_GetIndElem_Predicate(aIndElem,
          dfm2::CIsBV_IntersectLine<dfm2::CBV3_Sphere<double>>(src1,dir1),
          0, aNodeBVH, aAABB);
      delfem2::CPointElemSurf<double> pes;
      if( !aIndElem.empty() ) {
        std::map<double, delfem2::CPointElemSurf<double>> mapDepthPES;
        IntersectionRay_MeshTri3DPart(
            mapDepthPES,
            delfem2::CVec3d(src1),
            delfem2::CVec3d(dir1),
            aTri, aXYZ, aIndElem, 1.0e-10);
        if( !mapDepthPES.empty() ) {
          pes = mapDepthPES.begin()->second;
        }
      }
      if( pes.itri == UINT_MAX ) {
        aRGB[(ih * nwidth + iw) * 3 + 0] = 200;
        aRGB[(ih * nwidth + iw) * 3 + 1] = 255;
        aRGB[(ih * nwidth + iw) * 3 + 2] = 255;
      }
      else {
        const unsigned int itri = pes.itri;
        assert( itri < aTri.size()/3 );
        double n[3], area; delfem2::UnitNormalAreaTri3(
            n,area,
            aXYZ.data()+aTri[itri*3+0]*3,
            aXYZ.data()+aTri[itri*3+1]*3,
            aXYZ.data()+aTri[itri*3+2]*3);
        const double dot = n[0]*dir1[0] + n[1]*dir1[1] + n[2]*dir1[2];
        aRGB[(ih * nwidth + iw) * 3 + 0] = static_cast<unsigned char>(-dot*255);
        aRGB[(ih * nwidth + iw) * 3 + 1] = static_cast<unsigned char>(-dot*255);
        aRGB[(ih * nwidth + iw) * 3 + 2] = static_cast<unsigned char>(-dot*255);
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
    double rad = dfm2::CentsMaxRad_MeshTri3(aCent,
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
    dfm2::CLeafVolumeMaker_Mesh<dfm2::CBV3_Sphere<double>,double> lvm(
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
    tex.w = 200;
    tex.h = 200;
    tex.aRGB.resize(tex.w*tex.h*3);
    tex.max_x = +1.0;
    tex.min_x = -1.0;
    tex.max_y = +1.0;
    tex.min_y = -1.0;
    tex.z = 0.0;
  }

  dfm2::opengl::CViewer_GLFW viewer;
  viewer.width = 400;
  viewer.height = 400;
  viewer.nav.camera.view_height = 2;
  viewer.nav.camera.camera_rot_mode = dfm2::CCamera<double>::CAMERA_ROT_MODE::TBALL;

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
      float mMV[16], mP[16];
      viewer.nav.Mat4_MVP_OpenGL(mMV,mP,viewer.window);
      RayTracing(tex.aRGB,
           tex.h,tex.w, mMV,
           aNodeBVH,aAABB,aXYZ,aTri);
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


