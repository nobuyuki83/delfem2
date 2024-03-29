/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <iostream>
#include <vector>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/srch_trimesh3_class.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/isrf_adf.h"
#include "delfem2/msh_normal.h"
#include "delfem2/msh_unindexed.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_primitive.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/adf3.h"

namespace dfm2 = delfem2;

// -----------------------

// ---------------

void Draw(
    const delfem2::AdaptiveDistanceField3 &adf,
    const std::vector<double> &tri_xyz,
    const std::vector<double> &edge_xyz,
    bool is_show_cage) {
  const double color_face[3] = {1,1,1};
  //    std::cout << "ADF" << aNode.size() << std::endl;
  const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  if (!adf.aNode.empty() && is_show_cage) {
    delfem2::opengl::DrawThisAndChild_Wire(adf.aNode[0], adf.aNode);
  }

  if (!tri_xyz.empty()) {
    ::glColor3dv(color_face);
    ::glEnableClientState(GL_VERTEX_ARRAY);
    ::glVertexPointer(3, GL_DOUBLE, 0, tri_xyz.data() );
    ::glDrawArrays(GL_TRIANGLES, 0, tri_xyz.size()/3);
    ::glDisableClientState(GL_VERTEX_ARRAY);
  }
  if (!edge_xyz.empty() ){
    ::glColor3d(0, 0, 0);
    ::glEnableClientState(GL_VERTEX_ARRAY);
    ::glVertexPointer(3, GL_DOUBLE, 0, edge_xyz.data() );
    ::glDrawArrays(GL_LINES, 0, edge_xyz.size() / 3);
    ::glDisableClientState(GL_VERTEX_ARRAY);
  }
  if (is_lighting) { ::glEnable(GL_LIGHTING); }
}


std::vector<unsigned int> aTri;
std::vector<double> aXYZ;

void SetProblem(
    delfem2::AdaptiveDistanceField3 &adf,
    std::vector<double> &tri_xyz,
    std::vector<double> &edge_xyz,
    int iprob) {
  if (iprob == 0) {
    class CInSphere : public delfem2::Input_AdaptiveDistanceField3 {
     public:
      [[nodiscard]] double sdf(double x, double y, double z) const override {
        double n[3];
        return obj.Projection(n,
                              x, y, z);
      }
     public:
      delfem2::CSphere<double> obj;
    };
    CInSphere sphere;
    sphere.obj.radius_ = 0.5;
//    sphere.sdf.GetMesh(aTri, aXYZ, 0.01);
    double bb[6] = {-1, 1, -1, 1, -1, 1};
    adf.SetUp(sphere, bb);
    adf.BuildIsoSurface_MarchingCube(tri_xyz);
    delfem2::UnindexedEdgeMesh3_UnindexedTrimesh3(edge_xyz, tri_xyz);
  } else if (iprob == 1) {
    class CInTorus : public delfem2::Input_AdaptiveDistanceField3 {
     public:
      [[nodiscard]] double sdf(double x, double y, double z) const override {
        double n[3];
        return obj.Projection(n,
                              x, y, z);
      }
     public:
      delfem2::CTorus<double> obj;
    };
    CInTorus torus;
    torus.obj.radius_ = 0.5;
    torus.obj.radius_tube_ = 0.2;
//    torus.sdf.GetMesh(aTri, aXYZ, 0.01);
    double bb[6] = {-1, 1, -1, 1, -1, 1};
    adf.SetUp(torus, bb);
    adf.BuildIsoSurface_MarchingCube(tri_xyz);
    delfem2::UnindexedEdgeMesh3_UnindexedTrimesh3(edge_xyz, tri_xyz);
  } else if (iprob == 2) {
    class CMesh : public delfem2::Input_AdaptiveDistanceField3 {
     public:
      [[nodiscard]] double sdf(double x, double y, double z) const override {
        dfm2::CVec3d n0;
        double sdf0 = obj.SignedDistanceFunction(
            n0,
            dfm2::CVec3d(x, y, z),
            aXYZ, aTri, aNorm);
        return sdf0;
      }
     public:
      std::vector<double> aNorm;
      dfm2::CBVH_MeshTri3D<dfm2::CBV3d_Sphere, double> obj;
    };
    CMesh mesh;
    {
      std::cout << PATH_INPUT_DIR << std::endl;
      delfem2::Read_Ply(
          aXYZ, aTri,
          std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.ply");
          delfem2::Normalize_Points3(aXYZ, 1.7);
      mesh.obj.Init(
          aXYZ.data(), aXYZ.size() / 3,
          aTri.data(), aTri.size() / 3,
          0.0);
      mesh.aNorm.resize(aXYZ.size());
      delfem2::Normal_MeshTri3D(
          mesh.aNorm.data(),
          aXYZ.data(), aXYZ.size() / 3,
          aTri.data(), aTri.size() / 3);
    }
    double bb[6] = {-1, 1, -1, 1, -1, 1};
    adf.SetUp(mesh, bb);
    adf.BuildIsoSurface_MarchingCube(tri_xyz);
    delfem2::UnindexedEdgeMesh3_UnindexedTrimesh3(edge_xyz, tri_xyz);
  }
}

// ------------------------------------------------

int main() {
  delfem2::AdaptiveDistanceField3 adf;
  std::vector<double> tri_xyz;
  std::vector<double> edge_xyz;
  bool is_show_cage = false;
  //
  delfem2::glfw::CViewer3 viewer(2.0);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  double time_last_update = -3;
  int iproblem = 0;
  while (!glfwWindowShouldClose(viewer.window)) {
    {
      const double time = glfwGetTime();
      if (time - time_last_update > 2) {
        SetProblem(
            adf, tri_xyz, edge_xyz,
            iproblem);
        iproblem = (iproblem + 1) % 3;
        time_last_update = glfwGetTime();
      }
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    if (iproblem == 0) {
        is_show_cage = false;
      Draw(adf,tri_xyz,edge_xyz,is_show_cage);
    } else if (iproblem == 1) {
        is_show_cage = true;
      Draw(adf,tri_xyz,edge_xyz,is_show_cage);
    } else if (iproblem == 2) {
//      opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
      Draw(adf,tri_xyz,edge_xyz,is_show_cage);
    }
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
