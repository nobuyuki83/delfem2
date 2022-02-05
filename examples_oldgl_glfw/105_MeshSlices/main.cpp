/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <set>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshuni.h"
#include "delfem2/vec3.h"
#include "delfem2/slice.h"
#include "delfem2/geo_tri.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ---------------------------

void myGlutDisplay(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<delfem2::CSliceTriMesh> &aCS,
    const std::vector<std::set<unsigned int> > &ReebGraphCS,
    const std::vector<dfm2::CVec3d> &aCG_CS ) {
  ::glEnable(GL_LIGHTING);
  delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri);

  ::glDisable(GL_LIGHTING);
  ::glColor3d(1, 0, 0);
  ::glLineWidth(5);
  for (auto &cs : aCS) {
    ::glBegin(GL_LINE_LOOP);
    for (const auto &seg : cs.aTriInfo) {
      double pA[3], pB[3];
      seg.Pos3D(pA, pB,
                aXYZ, aTri);
      ::glVertex3d(pA[0], pA[1], pA[2]);
    }
    ::glEnd();
  }

  ::glDisable(GL_DEPTH_TEST);

  for (size_t ics = 0; ics < ReebGraphCS.size(); ++ics) {
    ::glColor3d(0, 0, 0);
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    ::glVertex3d(aCG_CS[ics].x, aCG_CS[ics].y, aCG_CS[ics].z);
    ::glEnd();
  }
  for (size_t ics = 0; ics < ReebGraphCS.size(); ++ics) {
    for (auto itr = ReebGraphCS[ics].begin(); itr != ReebGraphCS[ics].end(); ++itr) {
      const unsigned int jcs = *itr;
      assert(jcs < aCS.size());
      assert(abs(aCS[ics].IndHeight() - aCS[jcs].IndHeight()) == 1);
      ::glColor3d(0, 0, 0);
      ::glLineWidth(3);
      ::glBegin(GL_LINES);
      ::glVertex3d(aCG_CS[ics].x, aCG_CS[ics].y, aCG_CS[ics].z);
      ::glVertex3d(aCG_CS[jcs].x, aCG_CS[jcs].y, aCG_CS[jcs].z);
      ::glEnd();
    }
  }
  ::glEnable(GL_DEPTH_TEST);
}

void Hoge(
    std::vector<double> &vtx_xyz,
    std::vector<unsigned int> &tri_adjtri,
    std::vector<delfem2::CSliceTriMesh> &aCS,
    std::vector<std::set<unsigned int> > &ReebGraphCS,
    std::vector<dfm2::CVec3d> &aCG_CS) {
  delfem2::Read_Ply(
      vtx_xyz, tri_adjtri,
      std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.ply");
  delfem2::Normalize_Points3(vtx_xyz);
  std::vector<unsigned int> aTriSuTri;
  ElSuEl_MeshElem(
      aTriSuTri,
      tri_adjtri.data(), tri_adjtri.size() / 3, delfem2::MESHELEM_TRI,
      vtx_xyz.size() / 3);

  std::vector<double> aHeight;
  aHeight.push_back(-0.3);
  aHeight.push_back(-0.2);
  aHeight.push_back(-0.1);
  aHeight.push_back(-0.0);
  aHeight.push_back(+0.1);
  aHeight.push_back(+0.2);
  aHeight.push_back(+0.3);
  const double nrm[3] = {0, 1, 0};
  const double org[3] = {0, 0, 0};
  {
    std::vector<double> aHeightVtx(vtx_xyz.size() / 3);
    for (size_t ip = 0; ip < vtx_xyz.size() / 3; ++ip) {
      double x0 = vtx_xyz[ip * 3 + 0] - org[0];
      double y0 = vtx_xyz[ip * 3 + 1] - org[1];
      double z0 = vtx_xyz[ip * 3 + 2] - org[2];
      aHeightVtx[ip] = nrm[0] * x0 + nrm[1] * y0 + nrm[2] * z0;
    }
    Slice_MeshTri3D_Heights(aCS,
                            aHeight,
                            aHeightVtx,
                            tri_adjtri, aTriSuTri);
  }
  MakeReebGraph(ReebGraphCS,
                aCS, tri_adjtri, aTriSuTri);
  assert(aCS.size() == ReebGraphCS.size());

  aCG_CS.resize(aCS.size());
  for (size_t ics = 0; ics < aCS.size(); ++ics) {
    const double h0 = aHeight[aCS[ics].IndHeight()];
    const double po[3] = {org[0] + nrm[0] * h0, org[1] + nrm[1] * h0, org[2] + nrm[2] * h0};
    double sum_area = 0.0;
    dfm2::CVec3d cg(0, 0, 0);
    for (auto &iseg : aCS[ics].aTriInfo) {
      double pA[3], pB[3];
      iseg.Pos3D(pA, pB,
                 vtx_xyz, tri_adjtri);
      double n0[3];
      dfm2::Normal_Tri3(
          n0,
          pA, pB, po);
      const double area0 = n0[0] * nrm[0] + n0[1] * nrm[1] + n0[2] * nrm[2];
      sum_area += area0;
      cg.p[0] += area0 * (po[0] + pA[0] + pB[0]) / 3.0;
      cg.p[1] += area0 * (po[1] + pA[1] + pB[1]) / 3.0;
      cg.p[2] += area0 * (po[2] + pA[2] + pB[2]) / 3.0;
    }
    cg /= sum_area;
    aCG_CS[ics] = cg;
  }
}

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  std::vector<delfem2::CSliceTriMesh> aCS;
  std::vector<std::set<unsigned int> > ReebGraphCS;
  std::vector<dfm2::CVec3d> aCG_CS;
  //
  delfem2::glfw::CViewer3 viewer(0.5);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  delfem2::opengl::setSomeLighting();

  Hoge(
      vtx_xyz, tri_vtx,
      aCS, ReebGraphCS, aCG_CS );

  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.DrawBegin_oldGL();
    myGlutDisplay(
        vtx_xyz, tri_vtx,
        aCS, ReebGraphCS, aCG_CS );
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
