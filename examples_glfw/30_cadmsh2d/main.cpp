/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <math.h>
#include "delfem2/cad2d.h"
#include "delfem2/vec3.h"

// ----

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/opengl/gl24_funcs.h"
#include "delfem2/opengl/gl4_funcs.h"
#include "delfem2/opengl/gl4_v23dtricad.h"
#include "delfem2/opengl/glfw_cam.h"

// -----------------------------------

CNav3D_GLFW nav;
CCad2D cad;
CMeshDynTri2D dmsh;
std::vector<double> aW;

CShader_Cad2D shdr_cad;
CShader_MeshDTri2D shdr_dmsh;

// -----------------------------------

void draw(GLFWwindow* window)
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  float mMV[16], mP[16]; nav.Matrix_MVP(mMV, mP, window);
  shdr_cad.Draw(mP, mMV, cad);
  shdr_dmsh.Draw(mP, mMV);
  
  
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

void callback_mouse_button(GLFWwindow* window, int button, int action, int mods)
{
  nav.Mouse(window,button,action,mods);
  float px, py; nav.PosMouse2D(px, py, window);
  cad.Pick(px, py, nav.camera.view_height);
}

void callback_cursor_position(GLFWwindow* window, double xpos, double ypos)
{
  nav.Motion(window,xpos,ypos);
  if( nav.ibutton == 0 ){
    float px0,py0, px1,py1; nav.PosMove2D(px0,py0, px1,py1, window);
    cad.DragPicked(px1,py1, px0,py0);
    shdr_cad.MakeBuffer(cad);
    // --
    std::vector<double> aXYVtx = cad.XY_VtxCtrl_Face(0);
    int nv = aXYVtx.size()/2;
    int np = dmsh.aVec2.size();
    for(int ip=0;ip<np;++ip){
      dmsh.aVec2[ip].x = 0.0;
      dmsh.aVec2[ip].y = 0.0;
      for(int iv=0;iv<nv;++iv){
        dmsh.aVec2[ip].x += aW[ip*nv+iv]*aXYVtx[iv*2+0];
        dmsh.aVec2[ip].y += aW[ip*nv+iv]*aXYVtx[iv*2+1];
      }
    }
    shdr_dmsh.MakeBuffer(dmsh.aVec2, dmsh.aETri);
  }
}

void callback_scroll(GLFWwindow* window, double xoffset, double yoffset)
{
  nav.camera.scale *= pow(1.01,yoffset);
}


int main(void)
{
  GLFWwindow* window = myGLFW_OpenWindow(800,600);
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, callback_resize);
  glfwSetKeyCallback(            window, callback_key);
  glfwSetMouseButtonCallback(    window, callback_mouse_button);
  glfwSetCursorPosCallback(      window, callback_cursor_position);
  glfwSetScrollCallback(         window, callback_scroll);
    
  // glad: load all OpenGL function pointers
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  shdr_cad.Compile();
  shdr_dmsh.Compile();
  {
    {
      std::vector<double> aXY = {-1,-1, +1,-1, +1,+1, -1,+1};
      cad.AddPolygon(aXY);
    }
    shdr_cad.MakeBuffer(cad);
    {
      std::vector<int> aFlgPnt, aFlgTri;
      CMesher_Cad2D mesher;
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
        MeanValueCoordinate2D(aW.data()+nv*ixy,
                              dmsh.aVec2[ixy].x, dmsh.aVec2[ixy].y,
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
  
  nav.camera.view_height = 1.5;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
    
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  while (!glfwWindowShouldClose(window)) { draw(window); }
#endif
  
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

