/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <math.h>
#include "delfem2/vec3.h"
//
#include "delfem2/cad2_dtri2.h"

// ----

#if defined(_MSC_VER)
#  include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#endif

#include "delfem2/opengl/gl_funcs.h"
#include "delfem2/opengl/glnew_funcs.h"
#include "delfem2/opengl/glnew_v23dtricad.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -------------------------------------

class CCadDtri_Viewer : public delfem2::opengl::CViewer_GLFW {
public:
  CCadDtri_Viewer(){
    {
      std::vector<double> aXY = {-1,-1, +1,-1, +1,+1, -1,+1};
      cad.AddPolygon(aXY);
    }
  }
  void InitGL(){
    shdr_cad.Compile();
    shdr_dmsh.Compile();
    {
      shdr_cad.MakeBuffer(cad);
      {
        std::vector<int> aFlgPnt, aFlgTri;
        delfem2::CMesher_Cad2D mesher;
        mesher.edge_length = 0.08;
        mesher.Meshing(dmsh, cad);
        shdr_dmsh.MakeBuffer(dmsh.aVec2, dmsh.aETri);
      }
      {
        std::vector<double> aXYVtx = cad.XY_VtxCtrl_Face(0);
        const int nxy = dmsh.aVec2.size();
        const int nv = aXYVtx.size()/2;
        aW.resize(nxy*nv);
        for(int ixy=0;ixy<nxy;++ixy){
          dfm2::MeanValueCoordinate2D(aW.data()+nv*ixy,
                                      dmsh.aVec2[ixy].x(), dmsh.aVec2[ixy].y(),
                                      aXYVtx.data(), aXYVtx.size()/2);
          double sum = 0.0;
          for(int iv=0;iv<nv;++iv){
            sum += aW[nv*ixy+iv];
          }
          assert( fabs(sum-1)<1.0e-10 );
        }
      }
    }
    
    shdr_cad.is_show_face = false;
    
    camera.view_height = 1.5;
    camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  }
  void mouse_press(const float src[3], const float dir[3]) override {
    cad.Pick(src[0], src[1], camera.view_height);
  }
  void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) override {
    if( nav.ibutton == 0 ){
      cad.DragPicked(src1[0],src1[1], src0[0],src0[1]);
      shdr_cad.MakeBuffer(cad);
        // --
      std::vector<double> aXYVtx = cad.XY_VtxCtrl_Face(0);
      unsigned int nv = aXYVtx.size()/2;
      int np = dmsh.aVec2.size();
      for(int ip=0;ip<np;++ip){
        dmsh.aVec2[ip].p[0] = 0.0;
        dmsh.aVec2[ip].p[1] = 0.0;
        for(int iv=0;iv<nv;++iv){
          dmsh.aVec2[ip].p[0] += aW[ip*nv+iv]*aXYVtx[iv*2+0];
          dmsh.aVec2[ip].p[1] += aW[ip*nv+iv]*aXYVtx[iv*2+1];
        }
      }
      shdr_dmsh.MakeBuffer(dmsh.aVec2, dmsh.aETri);
    }
  }
public:
  delfem2::CCad2D cad;
  delfem2::CMeshDynTri2D dmsh;
  std::vector<double> aW;
  
  delfem2::opengl::CShader_Cad2D shdr_cad;
  delfem2::opengl::CShader_MeshDTri2D shdr_dmsh;
};

CCadDtri_Viewer viewer;

// -----------------------------------

void draw(GLFWwindow* window)
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  int nw, nh; glfwGetFramebufferSize(window, &nw, &nh);
  const float asp = (float)nw/nh;
  float mMV[16], mP[16]; viewer.camera.Mat4_MVP_OpenGL(mMV, mP, asp);
  viewer.shdr_cad.Draw(mP, mMV, viewer.cad);
  viewer.shdr_dmsh.Draw(mP, mMV);
  
  glfwSwapBuffers(window);
  glfwPollEvents();
}

void callback_key(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
    glfwSetWindowShouldClose(window, GL_TRUE);
  }
}

void callback_resize(GLFWwindow* window, int width, int height)
{
  glViewport(0, 0, width, height);
}

int main()
{
  viewer.Init_newGL();
    
  // glad: load all OpenGL function pointers
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
  
  viewer.InitGL();

    
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

