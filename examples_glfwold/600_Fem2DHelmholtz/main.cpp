/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include "delfem2/ilu_mats.h"
#include "delfem2/fem_emats.h"
#include "delfem2/cad2_dtri2.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mats.h"
#include "delfem2/vecxitrsol.h"

// ---------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ---------------------------------------------------

int main(int argc,char* argv[])
{
  std::vector<unsigned int> aTri1;
  std::vector<double> aXY1;
  std::vector< std::vector<unsigned int> > aaIP;
  int ipCenter = -1;
  {
    dfm2::CCad2D cad2d;
    {
      double xy[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
      std::vector<double> aXY(xy,xy+8);
      cad2d.AddPolygon(aXY);
    }
    cad2d.AddVtxFace(0.0, 0.0, 0);
    dfm2::CMesher_Cad2D mshr;
    mshr.edge_length = 0.05;
    dfm2::CMeshDynTri2D dmsh;
    mshr.Meshing(dmsh, cad2d);
    dfm2::MeshTri2D_Export(aXY1, aTri1,
                           dmsh.aVec2, dmsh.aETri);
    
    aaIP.resize(4);
    aaIP[0] = mshr.IndPoint_IndEdge(0, true, cad2d);
    aaIP[1] = mshr.IndPoint_IndEdge(1, true, cad2d);
    aaIP[2] = mshr.IndPoint_IndEdge(2, true, cad2d);
    aaIP[3] = mshr.IndPoint_IndEdge(3, true, cad2d);
    ipCenter = 4;
  }
  std::vector<int> aBCFlag; // boundary condition flag
  std::vector<std::complex<double> > aCVal;
  dfm2::CMatrixSparse<std::complex<double> > mat_A;
  dfm2::CPreconditionerILU<std::complex<double> > ilu_A;
  std::vector<std::complex<double> > vec_b;
  // ----------------------
  {
    const unsigned int np = aXY1.size()/2;
    aCVal.assign(np, std::complex<double>(0.0));
    aBCFlag.resize(np,0);
    aBCFlag[ipCenter] = 1;
    aCVal[ipCenter] = 1;
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aTri1.data(), aTri1.size()/3, 3, aXY1.size()/2);
    dfm2::JArray_Sort(psup_ind, psup);
    mat_A.Initialize(np, 1, true);
    mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
    ilu_A.Initialize_ILU0(mat_A);
  }
  {
    const unsigned int np = aXY1.size()/2;
    const unsigned int nDoF = np;
    const double wave_length = 0.4;
    mat_A.SetZero();
    vec_b.assign(nDoF, 0.0);
    dfm2::MergeLinSys_Helmholtz_MeshTri2D(mat_A,vec_b.data(),
                                          wave_length,
                                          aXY1.data(),aXY1.size()/2,
                                          aTri1.data(),aTri1.size()/3,
                                          aCVal.data());
    for(auto & ipl : aaIP){
      dfm2::MergeLinSys_SommerfeltRadiationBC_Polyline2D(mat_A,vec_b.data(),
                                                         wave_length,
                                                         aXY1.data(),aXY1.size()/2,
                                                         ipl.data(),ipl.size(),
                                                         aCVal.data());
    }
    mat_A.SetFixedBC(aBCFlag.data());
    dfm2::setRHS_Zero(vec_b, aBCFlag,0);
    std::vector<std::complex<double> > vec_x;
    ilu_A.SetValueILU(mat_A);
    ilu_A.DoILUDecomp();
    vec_x.resize(vec_b.size());
    /*
     std::vector<double> aConv = Solve_PBiCGStab(vec_b.data(), vec_x.data(),
     1.0e-4, 400, mat_A, ilu_A);
     */
    /*
     std::vector<double> aConv = Solve_BiCGSTAB_Complex(vec_b, vec_x,
     1.0e-4,400, mat_A);
     */
    std::vector<double> aConv = Solve_PCOCG(vec_b.data(), vec_x.data(),
                                            1.0e-4, 400, mat_A, ilu_A);
    std::cout << aConv.size() << " " << aConv[ aConv.size()-1 ] << std::endl;
    
    for(size_t ic=0;ic<aConv.size();++ic){
      std::cout << ic << " " << aConv[ic] << std::endl;
    }
    //  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
    //
    dfm2::XPlusAY(aCVal,
                  nDoF,aBCFlag, std::complex<double>(1.0),vec_x);
  }
  std::vector<double> aVal(aCVal.size(),0.1);
  for(size_t ip=0;ip<aVal.size();++ip){ aVal[ip] = aCVal[ip].real(); }
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.nav.camera.view_height = 1.5;
  viewer.Init_oldGL();
  while(!glfwWindowShouldClose(viewer.window)){
    {
      static int iframe = 0;
      double time_cur = iframe*0.01;
      std::complex<double> rot(cos(time_cur),sin(time_cur));
      for(size_t ip=0;ip<aVal.size();++ip){
        aVal[ip] = (rot*aCVal[ip]).real();
      }
      iframe++;
    }
      // ------------
    viewer.DrawBegin_oldGL();
    {
      dfm2::opengl::DrawMeshTri2D_Edge(aTri1,aXY1);
      ::glPointSize(2);
      ::glColor3d(0,0,0);
      dfm2::opengl::DrawPoints2d_Points(aXY1);
      
      std::vector< std::pair<double,dfm2::CColor> > colorMap;
      //  makeHeatMap_BlueGrayRed(colorMap, -0.2, +0.2);
      dfm2::ColorMap_BlueCyanGreenYellowRed(colorMap, -0.2f, +0.2f);
      dfm2::opengl::DrawMeshTri2D_ScalarP1(aXY1.data(),aXY1.size()/2,
                                           aTri1.data(),aTri1.size()/3,
                                           aVal.data(),1,colorMap);
    }
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
  return 0;
}
