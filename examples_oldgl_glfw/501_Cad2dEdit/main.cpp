/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/cad2dtriv2.h"
#include "delfem2/cad2_dtri2.h"
#include <GLFW/glfw3.h>
#include <cmath>
#include <stack>

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif
// -------------------------------------

int main(int argc,char* argv[])
{
  class CCAD2DViewer : public delfem2::glfw::CViewer3 {
  public:
    CCAD2DViewer(){
      const double poly[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
      cad.AddPolygon(std::vector<double>(poly,poly+8));
    }
    //
    void mouse_press(const float src[3], const float dir[3]) override{
      cad.Pick(src[0], src[1], this->camera.view_height);
    }
    void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) override{
      cad.DragPicked(src1[0], src1[1], src0[0], src0[1]);
    }
    //
    void Draw(){
      DrawBegin_oldGL();
      delfem2::opengl::Draw_CCad2D(cad);
      SwapBuffers();
    }
    delfem2::CCad2D cad;
  };
  // --------------------
  CCAD2DViewer viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.5;
  delfem2::opengl::setSomeLighting();
  // --------------------
  while(!glfwWindowShouldClose(viewer.window)){
    viewer.Draw();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

