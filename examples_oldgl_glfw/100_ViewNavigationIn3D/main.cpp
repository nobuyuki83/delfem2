/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief demo for view navigation in 3D
 * @details this demo is for showing CViewer_GLFW functionalities
 */

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#include <cstdlib>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

int main()
{
  delfem2::glfw::CViewer3 viewer;
  bool is_color = true;
  {
    // setting keyboard callback
    const std::function<void(int,int)> rr = [&is_color](int key, [[maybe_unused]] int mods) {
      if (key == GLFW_KEY_SPACE) { is_color = !is_color; }
    };
    viewer.keypress_callbacks.push_back(rr);
  }
  delfem2::glfw::InitGLOld();
  viewer.InitGL();

  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    if( is_color ){
      ::glClearColor(0,0,0,1);
      ::glClear(GL_COLOR_BUFFER_BIT );
    }
    glBegin(GL_TRIANGLES);
    glColor3f(1.f, 0.f, 0.f);
    glVertex3f(-0.6f, -0.4f, 0.f);
    glColor3f(0.f, 1.f, 0.f);
    glVertex3f(0.6f, -0.4f, 0.f);
    glColor3f(0.f, 0.f, 1.f);
    glVertex3f(0.f, 0.6f, 0.f);
    glEnd();
    viewer.SwapBuffers();

    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
