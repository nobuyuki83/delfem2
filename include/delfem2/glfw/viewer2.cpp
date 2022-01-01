/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer2.h"

#include <cstdlib>
#include <cassert>

#if defined(_WIN32)  // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
#  define GL_SILENCE_DEPRECATION  // remove
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

#include <GLFW/glfw3.h>

#include "delfem2/mat4.h"

// ---------------
namespace delfem2::viewer2 {

// static delfem2::glfw::CViewer2 *pViewer2 = nullptr;

static void glfw_callback_key(
    GLFWwindow *window,
    int key,
    [[maybe_unused]] int scancode,
    int action,
    int mods) {
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GL_TRUE);
  }
  auto pViewer2 = static_cast<delfem2::glfw::CViewer2 *>(glfwGetWindowUserPointer(window));
  if (action == GLFW_PRESS) { pViewer2->key_press(key, mods); }
  else if (action == GLFW_RELEASE) { pViewer2->key_release(key, mods); }
  else {}
}

static void glfw_callback_resize(
    [[maybe_unused]] GLFWwindow *window, int width, int height) {
  glViewport(0, 0, width, height);
}

static void glfw_callback_mouse_button(
    GLFWwindow *window,
    int button,
    int action,
    int mods) {
#ifdef  IMGUI_VERSION
  if( ImGui::GetIO().WantCaptureMouse ){ return; }
#endif
  auto pViewer2 = static_cast<delfem2::glfw::CViewer2 *>(glfwGetWindowUserPointer(window));
  assert(pViewer2 != nullptr);
  int width, height;
  glfwGetWindowSize(window, &width, &height);
  //const float asp = static_cast<float>(width) / static_cast<float>(height);
  {
    ::delfem2::CMouseInput &nav = pViewer2->nav;
    nav.imodifier = mods;
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    nav.mouse_x = (2.0 * x - width) / width;
    nav.mouse_y = (height - 2.0 * y) / height;
    if (action == GLFW_PRESS) {  // mouse down
      nav.ibutton = button;
      nav.mouse_x_down = nav.mouse_x;
      nav.mouse_y_down = nav.mouse_y;
    }
  }
  if (action == GLFW_PRESS) {
    const CMat4f mMV = pViewer2->GetModelViewMatrix();
    const CMat4f mP = pViewer2->GetProjectionMatrix();
    float src[3], dir[3];
    pViewer2->nav.MouseRay(
        src, dir,
        (mP * mMV).data());
    pViewer2->mouse_press(src);
  }
  if( action == GLFW_RELEASE ){
    pViewer2->mouse_release();
    ::delfem2::CMouseInput &nav = pViewer2->nav;
    nav.ibutton = -1;
  }
}

static void glfw_callback_cursor_position(
    GLFWwindow *window,
    double xpos,
    double ypos) {
#ifdef  IMGUI_VERSION
  if( ImGui::GetIO().WantCaptureMouse ){ return; }
#endif
  auto pViewer2 = static_cast<delfem2::glfw::CViewer2 *>(glfwGetWindowUserPointer(window));
  assert(pViewer2 != nullptr);
  int width, height;
  glfwGetWindowSize(window, &width, &height);
  const float asp = static_cast<float>(width) / static_cast<float>(height);
  { // update nav
    ::delfem2::CMouseInput &nav = pViewer2->nav;
    const double mov_end_x = (2.0 * xpos - width) / width;
    const double mov_end_y = (height - 2.0 * ypos) / height;
    nav.dx = mov_end_x - nav.mouse_x;
    nav.dy = mov_end_y - nav.mouse_y;
    nav.mouse_x = mov_end_x;
    nav.mouse_y = mov_end_y;
  }
  if (pViewer2->nav.ibutton == GLFW_MOUSE_BUTTON_LEFT && pViewer2->nav.imodifier == GLFW_MOD_SHIFT) {
    ::delfem2::CMouseInput &nav = pViewer2->nav;
    const float si = 1.f / pViewer2->scale;
    pViewer2->trans[0] += static_cast<float>(nav.dx) * asp * si;
    pViewer2->trans[1] += static_cast<float>(nav.dy) * si;
    return;
  }
  if( pViewer2->nav.ibutton != -1 ){
    const CMat4f mMV = pViewer2->GetModelViewMatrix();
    const CMat4f mP = pViewer2->GetProjectionMatrix();
    float src0[3], src1[3], dir0[3], dir1[3];
    pViewer2->nav.RayMouseMove(
        src0, src1, dir0, dir1,
        (mP * mMV).data());
    pViewer2->mouse_drag(src0, src1);
  }
}

static void glfw_callback_scroll(
    [[maybe_unused]] GLFWwindow *window,
    [[maybe_unused]] double xoffset,
    double yoffset) {
#ifdef  IMGUI_VERSION
  if( ImGui::GetIO().WantCaptureMouse ){ return; }
#endif
  auto pViewer2 = static_cast<delfem2::glfw::CViewer2 *>(glfwGetWindowUserPointer(window));
  assert(pViewer2 != nullptr);
  pViewer2->scale *= powf(1.01f, float(yoffset));
}

}  // namespace delfem2::viewer2

void delfem2::glfw::CViewer2::OpenWindow() {
  // glfw window creation
  // --------------------
  this->window = glfwCreateWindow(int(width),
                                  int(height),
                                  title.c_str(),
                                  nullptr,
                                  nullptr);
  if (this->window == nullptr) {
    glfwTerminate();
  }

  glfwMakeContextCurrent(
      this->window);
  glfwSetWindowUserPointer(
      this->window, this);
  glfwSetFramebufferSizeCallback(
      this->window, delfem2::viewer2::glfw_callback_resize);
  glfwSetKeyCallback(
      this->window, delfem2::viewer2::glfw_callback_key);
  glfwSetMouseButtonCallback(
      this->window, delfem2::viewer2::glfw_callback_mouse_button);
  glfwSetCursorPosCallback(
      this->window, delfem2::viewer2::glfw_callback_cursor_position);
  glfwSetScrollCallback(
      this->window, delfem2::viewer2::glfw_callback_scroll);
}

void delfem2::glfw::CViewer2::DrawBegin_oldGL() const {
  ::glfwMakeContextCurrent(window);
  //::glClearColor(0.8f, 1.0f, 1.0f, 1.0f);
  ::glClearColor(
      static_cast<float>(this->bgcolor[0]),
      static_cast<float>(this->bgcolor[1]),
      static_cast<float>(this->bgcolor[2]),
      static_cast<float>(this->bgcolor[3]));
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  // glnew will skip compilling following section
#ifdef GL_PROJECTION

  { // make sure that the stack is clear
    int n0;
    ::glGetIntegerv(GL_MODELVIEW_STACK_DEPTH, &n0);
    int n1;
    ::glGetIntegerv(GL_PROJECTION_STACK_DEPTH, &n1);
    assert(n0 == 1 && n1 == 1);
  }

  const CMat4f mMV = this->GetModelViewMatrix();
  const CMat4f mP = this->GetProjectionMatrix();
  {
    const CMat4f mZ = CMat4f::ScaleXYZ(1, 1, -1);
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    ::glMultMatrixf((mZ * mP).transpose().data());
  }
  {
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
    ::glMultMatrixf(mMV.transpose().data());
  }
  ::glEnable(GL_NORMALIZE);
#endif

}

void delfem2::glfw::CViewer2::SwapBuffers() const {
  glfwSwapBuffers(this->window);
}

void delfem2::glfw::CViewer2::ExitIfClosed() const {
  if (!glfwWindowShouldClose(this->window)) { return; }
  glfwDestroyWindow(this->window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

std::array<float,16> delfem2::glfw::CViewer2::GetProjectionMatrix() const {
  int w0, h0;
  glfwGetWindowSize(window, &w0, &h0);
  float asp = static_cast<float>(w0)/h0;
  std::array<float,16> m;
  const float si = view_height / scale;
  delfem2::Mat4_AffineProjectionOrtho(
      m.data(),
      -asp * si, +asp * si,
      -1.f * si, +1.f * si,
      -1.f, +1.f);
  return m;
}

std::array<float,16> delfem2::glfw::CViewer2::GetModelViewMatrix() const {
  std::array<float,16> m;
  delfem2::Mat4_Identity(m.data());
  m[0 * 4 + 3] += this->trans[0];
  m[1 * 4 + 3] += this->trans[1];
  return m;
}
