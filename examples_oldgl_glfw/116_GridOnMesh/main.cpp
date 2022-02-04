/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <vector>
#include <string>
#include <cstdlib>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/srchuni_v3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshuni.h"
#include "delfem2/vec3.h"
#include "delfem2/geo_vec3.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v3q.h"

#ifndef M_PI
#  define M_PI 3.141592
#endif

namespace dfm2 = delfem2;

// ---------------------------------------------------

std::vector<double> vtx_xyz;
std::vector<unsigned int> tri_vtx;
std::vector<double> vtx_norm;
std::vector<unsigned int> tri_adjtri;

double pos2d_org_corner[4][2] = {
    {-0.2, -0.4},
    {+0.2, -0.4},
    {+0.2, -0.0},
    {-0.2, -0.0},
};
std::vector<dfm2::PointOnSurfaceMesh<double>> aPES0;

const unsigned int ndiv = 8;
std::vector<dfm2::PointOnSurfaceMesh<double>> aPES1;

// ---------------------------------------------------

void InitializeProblem() {
  {
    dfm2::Read_Ply(
        vtx_xyz, tri_vtx,
        std::filesystem::path(PATH_INPUT_DIR) / "bunny_2k.ply");
    double cx, cy, cz, wx, wy, wz;
    dfm2::CenterWidth_Points3(
        cx, cy, cz, wx, wy, wz,
        vtx_xyz);
    dfm2::Translate_Points3(
        vtx_xyz,
        -cx, -cy, -cz);
    double wm = wx;
    wm = (wx > wm) ? wx : wm;
    wm = (wy > wm) ? wy : wm;
    wm = (wz > wm) ? wz : wm;
    dfm2::Scale_PointsX(vtx_xyz,
                        1.5 / wm);
    dfm2::Rotate_Points3(
        vtx_xyz,
        -M_PI * 0.5, 0.0, 0.0);
  }
  {
    std::vector<unsigned int> elsup_ind, elsup;
    dfm2::JArray_ElSuP_MeshElem(
        elsup_ind, elsup,
        tri_vtx.data(), tri_vtx.size() / 3, 3,
        vtx_xyz.size() / 3);
    dfm2::ElSuEl_MeshElem(
        tri_adjtri,
        tri_vtx.data(), tri_vtx.size() / 3, 3,
        elsup_ind, elsup,
        3, 2, dfm2::noelElemFace_Tri);
    assert(tri_adjtri.size() == tri_vtx.size());
  }
}

void UpdateProblem() {
  {
    aPES0.clear();
    dfm2::CVec3d dir(0, 0, -1);
    for (auto o2d : pos2d_org_corner) {
      const dfm2::CVec3d o3d(o2d[0], o2d[1], 10);
      std::map<double, dfm2::PointOnSurfaceMesh<double>> mapDepthPES;
      IntersectionRay_MeshTri3(
          mapDepthPES,
          o3d, dir, tri_vtx, vtx_xyz,
          0.0);
      if (!mapDepthPES.empty()) {
        aPES0.push_back(mapDepthPES.begin()->second);
      }
      /*
      const std::vector<CPointElemSurf> aPES = IntersectionLine_MeshTri3D(
          o3d, dir,
          aTri, aXYZ);
      std::map<double,CPointElemSurf> mapPES;
      for(auto pes : aPES){
        CVector3 p0 = pes.Pos_Tri(aXYZ,aTri);
        mapPES.insert( std::make_pair((p0-o3d)*dir, pes) );
      }
      if( mapPES.empty() ) { continue; }
      aPES0.push_back(mapPES.begin()->second);
       */
    }
  }
  {
    vtx_norm.resize(vtx_xyz.size());
    dfm2::Normal_MeshTri3D(
        vtx_norm.data(),
        vtx_xyz.data(), vtx_xyz.size() / 3,
        tri_vtx.data(), tri_vtx.size() / 3);
  }
  {
    aPES1.clear();
    dfm2::CVec3d p0 = aPES0[0].PositionOnMeshTri3(vtx_xyz, tri_vtx);
    dfm2::CVec3d p1 = aPES0[1].PositionOnMeshTri3(vtx_xyz, tri_vtx);
    dfm2::CVec3d p2 = aPES0[2].PositionOnMeshTri3(vtx_xyz, tri_vtx);
    dfm2::CVec3d p3 = aPES0[3].PositionOnMeshTri3(vtx_xyz, tri_vtx);
    dfm2::CVec3d n0 = aPES0[0].UnitNormalOnMeshTri3(vtx_norm, tri_vtx);
    dfm2::CVec3d n1 = aPES0[1].UnitNormalOnMeshTri3(vtx_norm, tri_vtx);
    dfm2::CVec3d n2 = aPES0[2].UnitNormalOnMeshTri3(vtx_norm, tri_vtx);
    dfm2::CVec3d n3 = aPES0[3].UnitNormalOnMeshTri3(vtx_norm, tri_vtx);
    for (unsigned int i = 0; i < ndiv + 1; ++i) {
      for (unsigned int j = 0; j < ndiv + 1; ++j) {
        double ri = (double) i / ndiv;
        double rj = (double) j / ndiv;
        double r0 = (1 - ri) * (1 - rj);
        double r1 = ri * (1 - rj);
        double r2 = ri * rj;
        double r3 = (1 - ri) * rj;
        const dfm2::CVec3d pA = r0 * p0 + r1 * p1 + r2 * p2 + r3 * p3;
        const dfm2::CVec3d nA = r0 * n0 + r1 * n1 + r2 * n2 + r3 * n3;
        const std::vector<dfm2::PointOnSurfaceMesh<double>> aPES = IntersectionLine_MeshTri3(
            pA, nA,
            tri_vtx, vtx_xyz,
            0.0);
        std::map<double, dfm2::PointOnSurfaceMesh<double>> mapPES;
        for (auto pes : aPES) {
          dfm2::CVec3d p_intersec = pes.PositionOnMeshTri3(vtx_xyz, tri_vtx);
          mapPES.insert(std::make_pair(Distance3(p_intersec, pA), pes));
        }
        if (mapPES.empty()) {
          aPES1.emplace_back();
        } else {
          aPES1.push_back(mapPES.begin()->second);
        }
      }
    }
    assert(aPES1.size() == (ndiv + 1) * (ndiv + 1));
  }

}

// -----------------------------------------------------

void myGlutDisplay() {
  ::glEnable(GL_LIGHTING);
  { // ball
    ::glDisable(GL_TEXTURE_2D);
    float gray[4] = {1.0f, 0.0f, 0.0f, 1.f};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    float shine[4] = {0, 0, 0, 0};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
    ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 127.0);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz, tri_vtx);
  }
  for (auto PES : aPES0) {
    auto v0 = PES.PositionOnMeshTri3(vtx_xyz, tri_vtx);
    float gray[4] = {0.0f, 0.0f, 1.0f, 1.f};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    dfm2::opengl::DrawSphereAt(32, 32, 0.02, v0[0], v0[1], v0[2]);
  }
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_DEPTH_TEST);
  ::glBegin(GL_QUADS);
  for (unsigned int idiv = 0; idiv < ndiv; ++idiv) {
    for (unsigned int jdiv = 0; jdiv < ndiv; ++jdiv) {
      dfm2::CVec3d p0 = aPES1[(idiv + 0) * (ndiv + 1) + (jdiv + 0)].PositionOnMeshTri3(vtx_xyz, tri_vtx);
      dfm2::CVec3d p1 = aPES1[(idiv + 1) * (ndiv + 1) + (jdiv + 0)].PositionOnMeshTri3(vtx_xyz, tri_vtx);
      dfm2::CVec3d p2 = aPES1[(idiv + 1) * (ndiv + 1) + (jdiv + 1)].PositionOnMeshTri3(vtx_xyz, tri_vtx);
      dfm2::CVec3d p3 = aPES1[(idiv + 0) * (ndiv + 1) + (jdiv + 1)].PositionOnMeshTri3(vtx_xyz, tri_vtx);
      if ((idiv + jdiv) % 2 == 0) {
        ::glColor3d(0, 0, 0);
      } else {
        ::glColor3d(1, 1, 1);
      }
      dfm2::opengl::myGlVertex(p0);
      dfm2::opengl::myGlVertex(p1);
      dfm2::opengl::myGlVertex(p2);
      dfm2::opengl::myGlVertex(p3);
    }
  }
  ::glEnd();
  ::glEnable(GL_DEPTH_TEST);
}

int main() {
  dfm2::glfw::CViewer3 viewer;
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  dfm2::opengl::setSomeLighting();

  InitializeProblem();
  UpdateProblem();

  int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window)) {
    iframe += 1;
    pos2d_org_corner[2][0] = +0.2 + 0.1 * sin(iframe * 0.2);
    pos2d_org_corner[2][1] = +0.0 + 0.1 * cos(iframe * 0.2);
    UpdateProblem();
    // -------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
