/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <vector>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/noise.h"

namespace dfm2 = delfem2;

// ---------------------

int main() {
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();

  unsigned int image_size = 256;
  std::vector<double> value_image;
  while (!glfwWindowShouldClose(viewer.window)) {
    {
      static int iframe = 0;
      if (iframe == 0) {
        dfm2::ComputePerlin(
            value_image,
            image_size, image_size,
            4, 4, 0.9);
      }
      iframe = (iframe + 1) % 30;
    }

    viewer.DrawBegin_oldGL();

    ::glColor3d(1, 1, 1);
    ::glBegin(GL_LINE_LOOP);
    ::glVertex2d(-0.5, -0.5);
    ::glVertex2d(+0.5, -0.5);
    ::glVertex2d(+0.5, +0.5);
    ::glVertex2d(-0.5, +0.5);
    ::glEnd();

    std::vector<std::pair<double, dfm2::CColor> > color_map;
    ColorMap_BlueCyanGreenYellowRed(color_map, -0.5, +0.5);
    //  makeHeatMap_BlueGrayRed(colorMap, -0.8, +0.8);
    ::glBegin(GL_QUADS);
    const unsigned int image_height = image_size;
    const unsigned int image_width = image_size;
    for (unsigned int jh = 0; jh < image_height - 1; ++jh) {
      for (unsigned int jw = 0; jw < image_width - 1; ++jw) {
        const unsigned int i00 = (jh + 0) * image_width + (jw + 0);
        const unsigned int i10 = (jh + 0) * image_width + (jw + 1);
        const unsigned int i11 = (jh + 1) * image_width + (jw + 1);
        const unsigned int i01 = (jh + 1) * image_width + (jw + 0);
        double v00 = value_image[i00];
        double v10 = value_image[i10];
        double v11 = value_image[i11];
        double v01 = value_image[i01];
        double x0 = -0.5 + 1.0 / (image_width - 1) * (jw + 0);
        double x1 = -0.5 + 1.0 / (image_width - 1) * (jw + 1);
        double y0 = -0.5 + 1.0 / (image_height - 1) * (jh + 0);
        double y1 = -0.5 + 1.0 / (image_height - 1) * (jh + 1);
        dfm2::opengl::heatmap(v00, color_map);
        ::glVertex2d(x0, y0);
        dfm2::opengl::heatmap(v10, color_map);
        ::glVertex2d(x1, y0);
        dfm2::opengl::heatmap(v11, color_map);
        ::glVertex2d(x1, y1);
        dfm2::opengl::heatmap(v01, color_map);
        ::glVertex2d(x0, y1);
      }
    }
    ::glEnd();

    ::glPointSize(5);

    // yerrow: input
    ::glLineWidth(1);
    ::glPointSize(5);
    ::glColor3d(1, 1, 0);

    // magenda: last
    ::glLineWidth(1);
    ::glPointSize(5);
    ::glColor3d(1, 0, 1);

    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
