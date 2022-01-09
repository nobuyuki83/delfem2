/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cassert>
#include <vector>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/vec3.h"
#include "delfem2/geo_curve_cubic.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

void Draw(
    const std::vector<dfm2::CVec3d> &cps,
    const std::vector<dfm2::CVec3d> &sample) {
  ::glPointSize(5);
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_POINTS);
  for (const auto &cp: cps) {
    ::glVertex3dv(cp.data());
  }
  ::glEnd();
  //
  ::glLineWidth(1);
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_LINE_STRIP);
  for (const auto &cp: cps) {
    ::glVertex3dv(cp.data());
  }
  ::glEnd();
  //
  ::glLineWidth(1);
  ::glColor3d(1, 0, 0);
  ::glBegin(GL_LINE_STRIP);
  for (const auto &p: sample) {
    ::glVertex3dv(p.data());
  }
  ::glEnd();
}

// Suppose that the knitted peice start from point(0,0) and stay in x-y plane
template<typename VEC>
std::vector<std::vector<VEC>> RepeatKnitUnit_RowIndependent(
    std::vector<std::vector<VEC>> &cps_knit,
    unsigned int row, unsigned int col,
    const std::vector<VEC> &cps_unit,
    const double stride_y) {
  double min_x = cps_unit[0][0], max_x = cps_unit[0][0];
  for (unsigned int i = 1; i < cps_unit.size(); i++) {
    min_x = min_x < cps_unit[i][0] ? min_x : cps_unit[i][0];
    max_x = max_x > cps_unit[i][0] ? min_x : cps_unit[i][0];
  }

  const double interval_x = max_x - min_x;
  const double interval_z = 0.;

  for (unsigned int j = 0; j < row; j++) {
    const double add_y = j * stride_y;
    const double add_z = j * interval_z;
    std::vector<VEC> yarn;
    yarn.reserve(cps_unit.size() * col);
    yarn.emplace_back(cps_unit[0][0], cps_unit[0][1] + add_y, cps_unit[0][2] + add_z);
    for (unsigned int i = 0; i < col; i++) {
      const double add_x = i * interval_x;
      for (unsigned int k = 1; k < cps_unit.size(); k++) { // ignore the first node to avoid repetation
        yarn.emplace_back(VEC(cps_unit[k][0] + add_x, cps_unit[k][1] + add_y, cps_unit[k][2] + add_z));
      }
    }
    cps_knit.emplace_back(yarn);
  }
  return cps_knit;
}

template<typename VEC>
void RepeatKnitUnit_SingleStrand(
    std::vector<VEC> &cps_knit,
    unsigned int row,
    unsigned int col,
    const std::vector<VEC> &cps,
    const double stride_y) {

  double min_x = cps[0][0], max_x = cps[0][0];
  for (unsigned int i = 1; i < cps.size(); i++) {
    min_x = min_x < cps[i][0] ? min_x : cps[i][0];
    max_x = max_x > cps[i][0] ? min_x : cps[i][0];
  }

  const double stride_x = max_x - min_x;
  const double interval_z = 0.;

  const unsigned int num_cp_unit = cps.size();
  for (unsigned int jrow = 0; jrow < row; jrow++) {
    const double add_y = jrow * stride_y;
    const double add_z = jrow * interval_z;
    std::vector<VEC> cps_row;
    cps_row.emplace_back(VEC(cps[0][0], cps[0][1] + add_y, cps[0][2] + add_z));
    for (unsigned int icol = 0; icol < col; icol++) {
      const double add_x = icol * stride_x;
      for (unsigned int k = 1; k < num_cp_unit; k++) { // ignore the first node to avoid repetition
        cps_row.emplace_back(VEC(cps[k][0] + add_x, cps[k][1] + add_y, cps[k][2] + add_z));
      }
    }
    if (jrow % 2 == 0) {
      for (unsigned int i = 0; i < cps_row.size(); i++) {
        cps_knit.emplace_back(cps_row.at(i));
      }
    } else {
      for (int i = static_cast<int>(cps_row.size() - 1); i >= 0; i--) {
        cps_knit.emplace_back(cps_row.at(i));
      }
    }
  }
}

int main() {
  std::vector<dfm2::CVec3d> cps_unit = {
      {-0.4, 0, 0},
      {-0.06, 0, 0},
      {-0.02, 0.45, 0.5},
      {-0.34, 0.45, 0.5},
      {-0.31, 1, 0.2},
      {0.31, 1, 0.2},
      {0.34, 0.45, 0.5},
      {0.02, 0.45, 0.5},
      {0.06, 0, 0},
      {0.4, 0, 0}};

  delfem2::glfw::CViewer3 viewer(1);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  {
    std::vector<dfm2::CVec3d> sample;
    delfem2::Polyline_CubicBSplineCurve(
        sample,
        1001, cps_unit);
    for (unsigned int iframe = 0; iframe < 100; ++iframe) {
      viewer.DrawBegin_oldGL();
      Draw(cps_unit, sample);
      ::glPushMatrix();
      ::glTranslated(0, 0.4, 0);
      Draw(cps_unit, sample);
      ::glPopMatrix();
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }
  { // RepeatKnitUnit_SingleStrand geometry. Each row separated
    std::vector<std::vector<delfem2::CVec3d>> cps_knit;
    RepeatKnitUnit_RowIndependent<delfem2::CVec3d>(
        cps_knit,
        3, 3, cps_unit, 0.5);
    std::vector<std::vector<dfm2::CVec3d >> samples;
    samples.resize(cps_knit.size());
    for (unsigned int iy = 0; iy < cps_knit.size(); ++iy) {
      delfem2::Polyline_CubicBSplineCurve(
          samples[iy],
          1001, cps_knit[iy]);
    }
    for (unsigned int iframe = 0; iframe < 100; iframe++) {
      viewer.DrawBegin_oldGL();
      for (unsigned int iy = 0; iy < cps_knit.size(); ++iy) {
        Draw(cps_knit[iy], samples[iy]);
      }
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }
  { // RepeatKnitUnit_SingleStrand geometry. Each row separated
    std::vector<delfem2::CVec3d> cps_knit;
    RepeatKnitUnit_SingleStrand<delfem2::CVec3d>(
        cps_knit,
        3, 3, cps_unit, 0.5);
    std::vector<dfm2::CVec3d> sample;
    delfem2::Polyline_CubicBSplineCurve(
        sample,
        1001, cps_knit);
    for (unsigned int iframe = 0; iframe < 100; iframe++) {
      viewer.DrawBegin_oldGL();
      Draw(cps_knit, sample);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
