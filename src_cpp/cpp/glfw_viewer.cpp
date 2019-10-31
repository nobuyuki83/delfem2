//
//  glfw_viewer.cpp
//  00_openwin
//
//  Created by Nobuyuki Umetani on 2019-10-31.
//

#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>

#include "delfem2/glfw_viewer.hpp"

// ---------------

static CViewer_GLFW* pViewer;

// ---------------

static void glfw_callback_error(int error, const char* description)
{
  fputs(description, stderr);
}

static void glfw_callback_key(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);
}

static void glfw_callback_resize(GLFWwindow* window, int width, int height)
{
  glViewport(0, 0, width, height);
}

static void glfw_callback_mouse_button(GLFWwindow* window, int button, int action, int mods)
{
  pViewer->nav.Mouse(window,button,action,mods);
}

static void glfw_callback_cursor_position(GLFWwindow* window, double xpos, double ypos)
{
  pViewer->nav.Motion(window,xpos,ypos);
}

static void glfw_callback_scroll(GLFWwindow* window, double xoffset, double yoffset)
{
  pViewer->nav.camera.scale *= pow(1.01,yoffset);
}


void CViewer_GLFW::Init_GLold()
{
  pViewer = this;
  // -----
  glfwSetErrorCallback(glfw_callback_error);
  if (!glfwInit()){
    exit(EXIT_FAILURE);
  }
  window = glfwCreateWindow(640, 480, "Simple example", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, glfw_callback_resize);
  glfwSetKeyCallback(            window, glfw_callback_key);
  glfwSetMouseButtonCallback(    window, glfw_callback_mouse_button);
  glfwSetCursorPosCallback(      window, glfw_callback_cursor_position);
  glfwSetScrollCallback(         window, glfw_callback_scroll);
}


void CViewer_GLFW::Init_GLnew()
{
  pViewer = this;
  // -----
  glfwSetErrorCallback(callback_error);
  if (!glfwInit()){
    exit(EXIT_FAILURE);
  }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  
#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif
  
  /*
   // Decide GL+GLSL versions
   #if __APPLE__
   // GL 3.2 + GLSL 150
   const char *glsl_version = "#version 150";
   glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
   glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
   glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // 3.2+ only
   glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);       // Required on Mac
   #else
   // GL 3.0 + GLSL 130
   const char *glsl_version = "#version 130";
   glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
   glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
   //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
   //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
   #endif
   */
  
    // glfw window creation
    // --------------------
  this->window = glfwCreateWindow(640,
                                  480,
                                  "LearnOpenGL",
                                  NULL,
                                  NULL);
  if (this->window == NULL)
  {
    glfwTerminate();
  }
  
  glfwMakeContextCurrent(this->window);
  glfwSetFramebufferSizeCallback(this->window, glfw_callback_resize);
  glfwSetKeyCallback(            this->window, glfw_callback_key);
  glfwSetMouseButtonCallback(    this->window, glfw_callback_mouse_button);
  glfwSetCursorPosCallback(      this->window, glfw_callback_cursor_position);
  glfwSetScrollCallback(         this->window, glfw_callback_scroll);
  return;
}


void CViewer_GLFW::DrawBegin_Glold()
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  float mP[16], mMV[16];
  nav.Matrix_MVP(mMV, mP, this->window);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMultMatrixf(mP);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixf(mMV);
}
