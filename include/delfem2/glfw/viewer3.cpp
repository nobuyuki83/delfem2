/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#if defined(_WIN32) // windows
  #define NOMINMAX   // to remove min,max macro
  #include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
  #define GL_SILENCE_DEPRECATION
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif

#include "delfem2/glfw/viewer3.h"
#include <cstdlib>
#include <cassert>
#include <GLFW/glfw3.h>

// ---------------

namespace delfem2 {
namespace glfw {
namespace viewer3 {

static delfem2::glfw::CViewer3 *pViewer3 = 0;

static void glfw_callback_key(GLFWwindow *window, int key, int scancode, int action, int mods) {
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GL_TRUE);
  }
  if (action == GLFW_PRESS) { pViewer3->key_press(key, mods); }
  else if (action == GLFW_RELEASE) { pViewer3->key_release(key, mods); }
}

static void glfw_callback_resize(GLFWwindow *window, int width, int height) {
  glViewport(0, 0, width, height);
}

static void glfw_callback_mouse_button(GLFWwindow *window, int button, int action, int mods) {
  assert(pViewer3 != 0);
  int width, height;
  glfwGetWindowSize(window, &width, &height);
  const float asp = width / (float) height;
  { // save input
    ::delfem2::CMouseInput &nav = pViewer3->nav;
    nav.imodifier = mods;
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    nav.mouse_x = (2.0 * x - width) / width;
    nav.mouse_y = (height - 2.0 * y) / height;
    if (action == GLFW_RELEASE) {
      nav.ibutton = -1;
    } else if (action == GLFW_PRESS) { // mouse down
      nav.ibutton = button;
      nav.mouse_x_down = nav.mouse_x;
      nav.mouse_y_down = nav.mouse_y;
    }
  }
  if (action == GLFW_PRESS && (mods == GLFW_MOD_SHIFT || mods == GLFW_MOD_ALT)) {
    // view control
    return;
  }
  if (action == GLFW_PRESS) { // "press callback"
    float src[3], dir[3];
    float mMVP[16];
    {
      float mMV[16], mP[16];
      pViewer3->camera.Mat4_MVP_OpenGL(mMV, mP, asp);
      ::delfem2::MatMat4(mMVP, mMV, mP);
    }
    pViewer3->nav.MouseRay(src, dir, asp, mMVP);
    pViewer3->mouse_press(src, dir);
  }
  if (action == GLFW_RELEASE) { // "release callback"
    pViewer3->mouse_release();
  }
}

static void glfw_callback_cursor_position(GLFWwindow *window, double xpos, double ypos) {
  assert(pViewer3 != 0);
  int width, height;
  glfwGetWindowSize(window, &width, &height);
  const float asp = width / (float) height;
  { // update nav
    ::delfem2::CMouseInput &nav = pViewer3->nav;
    const double mov_end_x = (2.0 * xpos - width) / width;
    const double mov_end_y = (height - 2.0 * ypos) / height;
    nav.dx = mov_end_x - nav.mouse_x;
    nav.dy = mov_end_y - nav.mouse_y;
    nav.mouse_x = mov_end_x;
    nav.mouse_y = mov_end_y;
  }
  if (pViewer3->nav.ibutton == GLFW_MOUSE_BUTTON_LEFT ) { // drag for view control
    ::delfem2::CMouseInput &nav = pViewer3->nav;
    if (nav.imodifier == GLFW_MOD_ALT ) {
      pViewer3->camera.Rot_Camera(nav.dx, nav.dy);
      return;
    }
    else if (nav.imodifier == GLFW_MOD_SHIFT ) {
      pViewer3->camera.Pan_Camera(nav.dx, nav.dy);
      return;
    }
  }
  // drag call back
  if (pViewer3->nav.ibutton == 0) {
    float src0[3], src1[3], dir0[3], dir1[3];
    float mMVP[16];
    {
      float mMV[16], mP[16];
      pViewer3->camera.Mat4_MVP_OpenGL(mMV, mP, asp);
      ::delfem2::MatMat4(mMVP, mMV, mP);
    }
    pViewer3->nav.RayMouseMove(src0, src1, dir0, dir1, asp, mMVP);
    pViewer3->mouse_drag(src0, src1, dir0);
  }
}

static void glfw_callback_scroll(GLFWwindow *window, double xoffset, double yoffset) {
  assert(pViewer3 != nullptr);
  pViewer3->camera.scale *= pow(1.01, yoffset);
  pViewer3->mouse_wheel(yoffset);
}

}
}
}

void delfem2::glfw::CViewer3::InitGL()
{
  delfem2::glfw::viewer3::pViewer3 = this;
    // glfw window creation
    // --------------------
  this->window = glfwCreateWindow(width,
                                  height,
                                  "LearnOpenGL",
                                  nullptr,
                                  nullptr);
  if (this->window == nullptr)
  {
    glfwTerminate();
  }
  
  glfwMakeContextCurrent(this->window);
  glfwSetFramebufferSizeCallback(this->window, delfem2::glfw::viewer3::glfw_callback_resize);
  glfwSetKeyCallback(            this->window, delfem2::glfw::viewer3::glfw_callback_key);
  glfwSetMouseButtonCallback(    this->window, delfem2::glfw::viewer3::glfw_callback_mouse_button);
  glfwSetCursorPosCallback(      this->window, delfem2::glfw::viewer3::glfw_callback_cursor_position);
  glfwSetScrollCallback(         this->window, delfem2::glfw::viewer3::glfw_callback_scroll);
}


void delfem2::glfw::CViewer3::DrawBegin_oldGL() const
{
  ::glfwMakeContextCurrent(window);
  //::glClearColor(0.8f, 1.0f, 1.0f, 1.0f);
  ::glClearColor(
      this->bgcolor[0],
      this->bgcolor[1],
      this->bgcolor[2],
      this->bgcolor[3]);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );

  // glnew will skip compilling following section
#ifdef GL_PROJECTION

  { // make sure that the stack is clear
    int n0; ::glGetIntegerv(GL_MODELVIEW_STACK_DEPTH,&n0);
    int n1; ::glGetIntegerv(GL_PROJECTION_STACK_DEPTH,&n1);
    assert( n0 == 1 && n1 == 1 );
  }

  float mMV[16], mP[16];
  {
    int width0, height0;
    glfwGetFramebufferSize(window, &width0, &height0);
    float asp = width0 / (float) height0;
    camera.Mat4_MVP_OpenGL(mMV,mP, asp);
  }

  ::glEnable(GL_NORMALIZE); // GL_NORMALIZE is not defiend on the modern OpenGLae
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMultMatrixf(mP);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixf(mMV);
#endif
  
}

void delfem2::glfw::CViewer3::SwapBuffers() const
{
  glfwSwapBuffers(this->window);
}

void delfem2::glfw::CViewer3::ExitIfClosed() const
{
  if ( !glfwWindowShouldClose(this->window) ) { return; }
  glfwDestroyWindow(this->window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}