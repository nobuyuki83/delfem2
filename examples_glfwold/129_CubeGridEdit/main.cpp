/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <climits>
#include "delfem2/vec3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/camera.h"
#include "delfem2/gridcube.h"
// -------
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/gridcube_glold.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// -----------------------------------------

int main(int argc,char* argv[])
{
  // -------------
  class CViewer_CubeGrid : public dfm2::opengl::CViewer_GLFW
  {
  public:
    CViewer_CubeGrid(){
      aCubeGrid.emplace_back(0,0,0 );
      org = dfm2::CVec3d(0,0,0);
    }
    void mouse_press(const float src[3], const float dir[3]) override {
      dfm2::CVec3d offsym(0,0,0);
      if( imode_sym == 2 ){ offsym.p[2] = -elen*0.5; }
      double src_pick0[3] = {src[0],src[1],src[2]};
      double dir_pick0[3] = {dir[0],dir[1],dir[2]};
      double offsym0[3];   offsym.CopyTo(offsym0);
      Pick_CubeGrid(icube_picked, iface_picked,
                    src_pick0,dir_pick0, elen, offsym0, aCubeGrid);
      if( edit_mode == EDIT_ADD ){
        int ivx1,ivy1,ivz1;
        Adj_CubeGrid(ivx1,ivy1,ivz1,
                     icube_picked,iface_picked,aCubeGrid);
        Add_CubeGrid(aCubeGrid,ivx1,ivy1,ivz1);
        if( imode_sym == 1 ){ Add_CubeGrid(aCubeGrid,ivx1,ivy1,-ivz1-1); }
        if( imode_sym == 2 ){ Add_CubeGrid(aCubeGrid,ivx1,ivy1,-ivz1); }
      }
      if( edit_mode == EDIT_DELETE ){
        int ivx1 = aCubeGrid[icube_picked].ivx;
        int ivy1 = aCubeGrid[icube_picked].ivy;
        int ivz1 = aCubeGrid[icube_picked].ivz;
        Del_CubeGrid(aCubeGrid,ivx1,ivy1,ivz1);
        if( imode_sym == 1 ){ Del_CubeGrid(aCubeGrid,ivx1,ivy1,-ivz1-1); }
        if( imode_sym == 2 ){ Del_CubeGrid(aCubeGrid,ivx1,ivy1,-ivz1); }
      }
    }
    void key_press(int key, int mods) override {
      if( key ==  GLFW_KEY_A ){
        glfwSetWindowTitle(this->window, "Edit Mode: Add");
        edit_mode = EDIT_ADD;
      }
      if( key ==  GLFW_KEY_D ){
        glfwSetWindowTitle(this->window, "Edit Mode: Delete");
        edit_mode = EDIT_DELETE;
      }
    }
    void Draw() const{
      dfm2::CVec3d offsym(0,0,0);
      if( imode_sym == 2 ){ offsym.p[2] = -elen*0.5; }
      for(unsigned int ic=0;ic<aCubeGrid.size();++ic){
        delfem2::opengl::Draw_CubeGrid(ic==icube_picked, iface_picked, elen, org+offsym, aCubeGrid[ic]);
      }
    }
  public:
    int imode_sym = 2;
    std::vector<dfm2::CCubeGrid> aCubeGrid;
    const double elen = 1.0;
    dfm2::CVec3d org;
    unsigned int icube_picked = -1;
    int iface_picked = -1;
    enum EDIT_MODE {
      EDIT_NONE,
      EDIT_ADD,
      EDIT_DELETE,
    };
    EDIT_MODE edit_mode = EDIT_ADD;
  } viewer;
  viewer.Init_oldGL();
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    viewer.Draw();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


