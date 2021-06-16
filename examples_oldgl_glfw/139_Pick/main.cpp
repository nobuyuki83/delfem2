/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <algorithm>

namespace dfm2 = delfem2;

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ; // 3d points
  std::vector<unsigned int> aTri;

  { // load input mesh
    const auto path =std::string(PATH_SOURCE_DIR) + "/../../test_inputs/arm_16k.ply";
    std::cout << path << std::endl;
    delfem2::Read_Ply(
        path,
        aXYZ, aTri);
    dfm2::Normalize_Points3(aXYZ, 2.0);
    std::cout << "point_size: " << aXYZ.size()/3 << std::endl;
    std::cout << "triangle_size: " << aTri.size()/3 << std::endl;
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

    dfm2::CLeafVolumeMaker_Mesh<dfm2::CBV3_Sphere<double>,double> lvm(
        1.0e-10,
        aXYZ.data(), aXYZ.size()/3,
        aTri.data(), aTri.size()/3, 3);
    dfm2::BVH_BuildBVHGeometry(aAABB,
        0, aNodeBVH,
        lvm);
    { // assertion
      dfm2::Check_MortonCode_Sort(aSortedId, aSortedMc, aCent, min_xyz, max_xyz);
      dfm2::Check_MortonCode_RangeSplit(aSortedMc);
      dfm2::Check_BVH(aNodeBVH,aCent.size()/3);
    }
  }

  std::vector<unsigned int> aFlagElem(aTri.size()/3, 0);

  class MyView
      : public dfm2::glfw::CViewer3
  {
  public:
    MyView(
        std::vector<unsigned int>& aFlagElem0,
        const std::vector<dfm2::CNodeBVH2>& aNodeBVH0,
        const std::vector<dfm2::CBV3_Sphere<double>>& aAABB0)
        : aFlagElem(aFlagElem0), aNodeBVH(aNodeBVH0), aAABB(aAABB0) {}
    void mouse_press(const float src[3], const float dir[3]) override {
      std::vector<unsigned int> aIndElem;
      dfm2::BVH_GetIndElem_Predicate(
          aIndElem,
          delfem2::CIsBV_IntersectLine< dfm2::CBV3_Sphere<double>, float>(src,dir),
          0, aNodeBVH, aAABB);
      for(auto ie : aIndElem) {
        aFlagElem[ie] = 1;
      }
    }
  public:
    std::vector<unsigned int>& aFlagElem;
    const std::vector<dfm2::CNodeBVH2>& aNodeBVH;
    const std::vector<dfm2::CBV3_Sphere<double>>& aAABB;
  };

  MyView viewer(
      aFlagElem,
      aNodeBVH,aAABB);

  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  const std::vector< std::pair<int,delfem2::CColor> > aColor = {
      { 2, delfem2::CColor::White() },
      { 2, delfem2::CColor::Red() } };

  while (!glfwWindowShouldClose(viewer.window))
  {
    // ----------
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    delfem2::opengl::DrawMeshTri3D_Edge(aXYZ,aTri);
    delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(
        aXYZ,aTri,aFlagElem,aColor);
    // ----------
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


