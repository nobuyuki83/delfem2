/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/lsilu_mats.h"
#include "delfem2/lsmats.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/femhelmholtz.h"
#include "delfem2/cad2.h"
#include "delfem2/cad2_mesher.h"
#include "delfem2/mshuni.h"
#include "delfem2/jagarray.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ---------------------------------------------------

int main() {
  std::vector<unsigned int> tri_vtx;
  std::vector<double> vtx_xy;
  std::vector<std::vector<unsigned int> > aaIP;
  int ipCenter;
  {
    dfm2::CCad2D cad2d;
    {
      double xy[8] = {-1, -1, +1, -1, +1, +1, -1, +1};
      std::vector<double> aXY(xy, xy + 8);
      cad2d.AddPolygon(aXY);
    }
    cad2d.AddVtxFace(0.0, 0.0, 0);
    dfm2::CMesher_Cad2D mshr;
    mshr.edge_length = 0.05;
    dfm2::CMeshDynTri2D dmsh;
    mshr.Meshing(dmsh, cad2d);
    dfm2::MeshTri2D_Export(
        vtx_xy, tri_vtx,
        dmsh.aVec2, dmsh.aETri);

    aaIP.resize(4);
    aaIP[0] = mshr.IndPoint_IndEdge(0, true, cad2d);
    aaIP[1] = mshr.IndPoint_IndEdge(1, true, cad2d);
    aaIP[2] = mshr.IndPoint_IndEdge(2, true, cad2d);
    aaIP[3] = mshr.IndPoint_IndEdge(3, true, cad2d);
    ipCenter = 4;
  }
  using COMPLEX = std::complex<double>;
  std::vector<int> aBCFlag; // boundary condition flag
  std::vector<COMPLEX> aCVal;
  dfm2::CMatrixSparse<COMPLEX> mat_A;
  dfm2::CPreconditionerILU<COMPLEX> ilu_A;
  std::vector<COMPLEX> vec_b;
  // ----------------------
  {
    const auto np = static_cast<unsigned int>(vtx_xy.size() / 2);
    aCVal.assign(np, COMPLEX(0.0));
    aBCFlag.resize(np, 0);
    aBCFlag[ipCenter] = 1;
    aCVal[ipCenter] = 1;
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(
        psup_ind, psup,
        tri_vtx.data(), tri_vtx.size() / 3, 3,
        vtx_xy.size() / 2);
    dfm2::JArray_Sort(psup_ind, psup);
    mat_A.Initialize(np, 1, true);
    mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
    ilu_A.SetPattern0(mat_A);
  }
  {
    const auto np = static_cast<unsigned int>(vtx_xy.size() / 2);
    const unsigned int nDoF = np;
    const double wave_length = 0.4;
    mat_A.setZero();
    vec_b.assign(nDoF, 0.0);
    dfm2::MergeLinSys_Helmholtz_MeshTri2D(
        mat_A, vec_b.data(),
        wave_length,
        vtx_xy.data(), vtx_xy.size() / 2,
        tri_vtx.data(), tri_vtx.size() / 3,
        aCVal.data());
    for (auto &ipl : aaIP) {
      dfm2::MergeLinSys_SommerfeltRadiationBC_Polyline2D(
          mat_A, vec_b.data(),
          wave_length,
          vtx_xy.data(), vtx_xy.size() / 2,
          ipl.data(), ipl.size(),
          aCVal.data());
    }
    mat_A.SetFixedBC(aBCFlag.data());
    dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
    std::vector<COMPLEX> vec_x;
    ilu_A.CopyValue(mat_A);
    ilu_A.Decompose();
    vec_x.resize(vec_b.size());
    /*
     std::vector<double> aConv = Solve_PBiCGStab(vec_b.data(), vec_x.data(),
     1.0e-4, 400, mat_A, ilu_A);
     */
    /*
     std::vector<double> aConv = Solve_BiCGSTAB_Complex(vec_b, vec_x,
     1.0e-4,400, mat_A);
     */
    std::vector<double> aConv = Solve_PCOCG(
        vec_b.data(), vec_x.data(),
        1.0e-4, 400, mat_A, ilu_A);
    std::cout << aConv.size() << " " << aConv[aConv.size() - 1] << std::endl;

    for (size_t ic = 0; ic < aConv.size(); ++ic) {
      std::cout << ic << " " << aConv[ic] << std::endl;
    }
    //  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
    //
    dfm2::XPlusAY(
        aCVal,
        nDoF, aBCFlag, std::complex<double>(1.0), vec_x);
  }
  std::vector<double> aVal(aCVal.size(), 0.1);
  for (size_t ip = 0; ip < aVal.size(); ++ip) { aVal[ip] = aCVal[ip].real(); }

  dfm2::glfw::CViewer3 viewer(1.5);
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  //
  while (!glfwWindowShouldClose(viewer.window)) {
    {
      static int iframe = 0;
      double time_cur = iframe * 0.01;
      std::complex<double> rot(cos(time_cur), sin(time_cur));
      for (size_t ip = 0; ip < aVal.size(); ++ip) {
        aVal[ip] = (rot * aCVal[ip]).real();
      }
      iframe++;
    }
    // ------------
    viewer.DrawBegin_oldGL();
    {
      dfm2::opengl::DrawMeshTri2D_Edge(tri_vtx, vtx_xy);
      ::glPointSize(2);
      ::glColor3d(0, 0, 0);
      dfm2::opengl::DrawPoints2d_Points(vtx_xy);

      std::vector<std::pair<double, dfm2::CColor> > colorMap;
      //  makeHeatMap_BlueGrayRed(colorMap, -0.2, +0.2);
      dfm2::ColorMap_BlueCyanGreenYellowRed(colorMap, -0.2f, +0.2f);
      dfm2::opengl::DrawMeshTri2D_ScalarP1(
          vtx_xy.data(), vtx_xy.size() / 2,
          tri_vtx.data(), tri_vtx.size() / 3,
          aVal.data(), 1, colorMap);
    }
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
