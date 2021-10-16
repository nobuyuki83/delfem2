/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>
#include <Eigen/LU>

#include "delfem2/vec2.h"
#include "delfem2/points.h"
#include "delfem2/mshuni.h"
#include "delfem2/slice.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/color.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"

template <class V2>
void solveSignFieldRBF (
    std::vector<double>& aWeight,
    std::vector<V2>& stroke_offset,
    const std::vector<V2>& stroke,
    double eps)
{
  const size_t n = stroke.size();
  std::vector<V2> aNorm;
  aNorm.resize(n, V2(0, 0));
  for(unsigned int is=0;is<stroke.size()-1;is++){
    V2 p0 = stroke[is+0];
    V2 p1 = stroke[is+1];
    V2 e = (p1-p0).normalized();
    V2 v(e.y, -e.x);
    aNorm[is+0] += v;
    aNorm[is+1] += v;
  }
  for(unsigned int i=0;i<n;i++){ aNorm[i].normalize(); }

  stroke_offset.resize(n*2);
  for(unsigned int i=0;i<n*2;i++){
    unsigned int i0 = i/2;
    V2 pi = stroke[i0];
    pi += (i%2==0) ? aNorm[i0]*eps : -1.0*aNorm[i0]*eps;
    stroke_offset[i] = pi;
  }
  // ---------
  Eigen::MatrixXd A(2*n,2*n);
  A.setZero();
  Eigen::VectorXd y(2*n);
  for(unsigned int i=0;i<n*2;i++){
    for(unsigned int j=0;j<n*2;j++){
      V2 pi = stroke_offset[i];
      V2 pj = stroke_offset[j];
      double r = (pi-pj).norm();
      A(i,j) = exp(-r);
//      if( i != j ) { A(i, j) = r * r * log(r); }
    }
    y(i) = (i%2==0) ? eps : -eps;
  }
  Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
  Eigen::VectorXd x = solver.solve(y);
  aWeight.resize(2*n);
  for(unsigned int i=0;i<2*n;i++){
    aWeight[i] = x(i);
  }
}

template <class V2>
double evaluateRBF(
    const V2& p,
    const std::vector<double>& aWeight,
    const std::vector<V2>& stroke_offset)
{
  int n = (int)aWeight.size();
  double t = 0;
  for(int i=0;i<n;i++){
    V2 pi = stroke_offset[i];
    double r = (p-pi).norm();
    t += aWeight[i]*exp(-r);
    //t += aWeight[i]*r*r*log(r);
  }
  return t;
}



// ---------------------------

void Draw(
    const std::vector<double> &vtx_xy,
    const std::vector<unsigned int> &tri_vtx,
    const std::vector<double> &vtx_val,
    const std::vector<delfem2::CSegInfo> &aSeg) {
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  delfem2::opengl::DrawMeshTri2D_Edge(
      vtx_xy.data(), vtx_xy.size() / 2,
      tri_vtx.data(), tri_vtx.size() / 3);
  std::vector<std::pair<double, delfem2::CColor> > colorMap;
  delfem2::ColorMap_BlueCyanGreenYellowRed(
      colorMap,
      -0.01, +0.01);
  delfem2::opengl::DrawMeshTri2D_ScalarP1(
      vtx_xy.data(), vtx_xy.size() / 2,
      tri_vtx.data(), tri_vtx.size() / 3,
      vtx_val.data(), 1,
      colorMap);
  ::glLineWidth(5);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (auto &iseg : aSeg) {
    double pA[2], pB[2];
    iseg.Pos2D(pA, pB,
               vtx_xy.data(), tri_vtx.data());
    ::glVertex2dv(pA);
    ::glVertex2dv(pB);
  }
  ::glEnd();

}

void Initialize(
    std::vector<double> &vtx_xy,
    std::vector<unsigned int> &tri_vtx,
    int ndiv) {
  std::vector<unsigned int> aQuad;
  delfem2::MeshQuad2D_Grid(
      vtx_xy, aQuad,
      ndiv, ndiv);
  delfem2::convert2Tri_Quad(
      tri_vtx,
      aQuad);
  delfem2::Translate_Points2(
      vtx_xy,
      -ndiv * 0.5, -ndiv * 0.5);
  delfem2::Scale_PointsX(
      vtx_xy,
      1.0 / ndiv);
}

int main() {
  std::vector<double> vtx_xy;
  std::vector<unsigned int> tri_vtx;
  Initialize(
      vtx_xy, tri_vtx,
      32);
  delfem2::glfw::CViewer3 viewer;

  //
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  while (!glfwWindowShouldClose(viewer.window)) {
    std::vector<delfem2::CSegInfo> segments;
    std::vector<double> vtx_val;
    {
      double time0 = glfwGetTime();
      auto ndiv = static_cast<unsigned int>(28. * sin(time0) + 30.0);
      std::vector<delfem2::CVec2d> stroke;
      for(unsigned int i=0;i<ndiv;++i){
        stroke.emplace_back(0.3*sin(i*0.1),0.3*cos(i*0.1));
      }
      std::vector<double> aWeight;
      std::vector<delfem2::CVec2d> stroke_offset;
      solveSignFieldRBF (
          aWeight,
          stroke_offset, stroke, 1.0e-3);
      vtx_val.resize(vtx_xy.size()/2);
      for(unsigned int ip=0;ip<vtx_xy.size()/2;++ip){
        vtx_val[ip] = evaluateRBF(
            delfem2::CVec2d(vtx_xy[ip*2+0], vtx_xy[ip*2+1]),
            aWeight, stroke_offset);
      }
      delfem2::AddContour(
          segments,
          0.0,
          tri_vtx.data(), tri_vtx.size() / 3,
          vtx_val.data());
    }
    // ------------------
    viewer.DrawBegin_oldGL();
    Draw(vtx_xy, tri_vtx, vtx_val, segments);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
