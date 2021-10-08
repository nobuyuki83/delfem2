/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer3.h"

#include <cstdlib>
#include <cassert>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

// ---------------

namespace delfem2::glfw::viewer3 {

static void glfw_callback_key(
    GLFWwindow *window,
    int key,
    [[maybe_unused]] int scancode,
    int action,
    int mods) {
  auto pViewer3 = static_cast<delfem2::glfw::CViewer3 *>(glfwGetWindowUserPointer(window));
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GL_TRUE);
  }
  if (action == GLFW_PRESS) {
    if (key == GLFW_KEY_PAGE_UP) { pViewer3->scale *= 1.03; }
    if (key == GLFW_KEY_PAGE_DOWN) { pViewer3->scale *= (1.0 / 1.03); }
    pViewer3->key_press(key, mods);
    for(const auto& func : pViewer3->keypress_callbacks){
      func(key, mods);
    }
  } else if (action == GLFW_RELEASE) { pViewer3->key_release(key, mods);
  } else if (action == GLFW_REPEAT) { pViewer3->key_repeat(key, mods); }

}

static void glfw_callback_resize(
    [[maybe_unused]] GLFWwindow *window,
    int width, int height) {
  glViewport(0, 0, width, height);
}

static void glfw_callback_mouse_button(
    GLFWwindow *window,
    int button,
    int action,
    int mods) {
  auto pViewer3 = static_cast<delfem2::glfw::CViewer3 *>(glfwGetWindowUserPointer(window));
//  std::cout << window << " key" << " " << pViewer3 << std::endl;
  assert(pViewer3 != nullptr);
  int width, height;
  glfwGetWindowSize(window, &width, &height);
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
    const CMat4f mP = pViewer3->GetProjectionMatrix();
    const CMat4f mMV = pViewer3->GetModelViewMatrix();
    {
      float src[3], dir[3];
      pViewer3->nav.MouseRay(
          src, dir,
          (mP * mMV).data());
      pViewer3->mouse_press(src, dir);
    }
  }
  if (action == GLFW_RELEASE) { // "release callback"
    pViewer3->mouse_release();
  }
}

static void glfw_callback_cursor_position(
    GLFWwindow *window,
    double xpos, double ypos) {
  auto pViewer3 = static_cast<delfem2::glfw::CViewer3 *>(glfwGetWindowUserPointer(window));
  assert(pViewer3 != nullptr);
  pViewer3->CursorPosition(xpos, ypos);
}

static void glfw_callback_scroll(
    GLFWwindow *window,
    [[maybe_unused]] double xoffset,
    double yoffset) {
  auto pViewer3 = static_cast<delfem2::glfw::CViewer3 *>(glfwGetWindowUserPointer(window));
  assert(pViewer3 != nullptr);
  pViewer3->scale *= pow(1.01, yoffset);
  pViewer3->mouse_wheel(yoffset);
}

}  // namespace delfem2

// ------------------------

std::array<float,16> delfem2::glfw::CViewer3::GetProjectionMatrix() const {
  int w0, h0;
  glfwGetWindowSize(window, &w0, &h0);
  const CMat4f mS = CMat4f::Scale((float)scale);
  const float asp = static_cast<float>(w0) / static_cast<float>(h0);
  const CMat4f mP = projection->GetMatrix(asp);
  return (mP * mS).GetStlArray();
}

void delfem2::glfw::CViewer3::CursorPosition(double xpos, double ypos) {
  int width0, height0;
  glfwGetWindowSize(window, &width0, &height0);
  { // update nav
    const double mov_end_x = (2.0 * xpos - width0) / width0;
    const double mov_end_y = (height0 - 2.0 * ypos) / height0;
    nav.dx = mov_end_x - nav.mouse_x;
    nav.dy = mov_end_y - nav.mouse_y;
    nav.mouse_x = mov_end_x;
    nav.mouse_y = mov_end_y;
  }
  if (this->nav.ibutton == GLFW_MOUSE_BUTTON_LEFT) {  // drag for view control
    if (nav.imodifier == GLFW_MOD_ALT) {
      this->view_rotation->Rot_Camera(
          static_cast<float>(nav.dx),
          static_cast<float>(nav.dy));
      return;
    } else if (nav.imodifier == GLFW_MOD_SHIFT) {
      const delfem2::CMat4f mP = this->GetProjectionMatrix();
      const float sx = (mP(3,3) - mP(0,3))/mP(0,0);
      const float sy = (mP(3,3) - mP(1,3))/mP(1,1);
      this->trans[0] += sx*nav.dx;
      this->trans[1] += sy*nav.dy;
      return;
    }
  }
  // drag call back
  if (this->nav.ibutton == 0) {
    const CMat4f mP = this->GetProjectionMatrix();
    const CMat4f mMV = this->GetModelViewMatrix();
    {
      float src0[3], src1[3], dir0[3], dir1[3];
      this->nav.RayMouseMove(
          src0, src1, dir0, dir1,
          (mP * mMV).data());
      this->mouse_drag(src0, src1, dir0);
    }
  }
}

void delfem2::glfw::CViewer3::InitGL() {
  namespace lcl = delfem2::glfw::viewer3;
  this->window = glfwCreateWindow(
      static_cast<int>(width),
      static_cast<int>(height),
      window_title.c_str(),
      nullptr,
      nullptr);
  if (this->window == nullptr) {
    glfwTerminate();
  }
  glfwMakeContextCurrent(this->window);
  glfwSetWindowUserPointer(this->window, this);
  glfwSetFramebufferSizeCallback(this->window, lcl::glfw_callback_resize);
  glfwSetKeyCallback(this->window, lcl::glfw_callback_key);
  glfwSetMouseButtonCallback(this->window, lcl::glfw_callback_mouse_button);
  glfwSetCursorPosCallback(this->window, lcl::glfw_callback_cursor_position);
  glfwSetScrollCallback(this->window, lcl::glfw_callback_scroll);
}

void delfem2::glfw::CViewer3::DrawBegin_oldGL() const {
  ::glfwMakeContextCurrent(window);
  //::glClearColor(0.8f, 1.0f, 1.0f, 1.0f);
  ::glClearColor(
      static_cast<GLclampf>(this->bgcolor[0]),
      static_cast<GLclampf>(this->bgcolor[1]),
      static_cast<GLclampf>(this->bgcolor[2]),
      static_cast<GLclampf>(this->bgcolor[3]));
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  // glnew will skip compilling following section
#ifdef GL_PROJECTION

#ifndef NDEBUG
  {  // make sure that the stack is clear
    int n0;
    ::glGetIntegerv(GL_MODELVIEW_STACK_DEPTH, &n0);
    int n1;
    ::glGetIntegerv(GL_PROJECTION_STACK_DEPTH, &n1);
    assert(n0 == 1 && n1 == 1);
  }
#endif

  ::glEnable(GL_NORMALIZE); // GL_NORMALIZE is not defiend on the modern OpenGL
  {
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    const CMat4f mP = this->GetProjectionMatrix();
    const CMat4f mZ = CMat4f::ScaleXYZ(1, 1, -1);
    ::glMultMatrixf((mZ * mP).transpose().data());
  }
  {
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
    const CMat4f mMV = this->GetModelViewMatrix();
    ::glMultMatrixf(mMV.transpose().data());
  }
#endif
}

void delfem2::glfw::CViewer3::SwapBuffers() const {
  glfwSwapBuffers(this->window);
}

void delfem2::glfw::CViewer3::ExitIfClosed() const {
  if (!glfwWindowShouldClose(this->window)) { return; }
  glfwDestroyWindow(this->window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
