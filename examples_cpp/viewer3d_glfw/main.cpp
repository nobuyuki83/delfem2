#include <iostream>
#include <math.h>

#if defined(_MSC_VER)
#include <windows.h> 
#endif

#include <GLFW/glfw3.h>

#include "delfem2/glfw_funcs.h"
#include "delfem2/gl_funcs.h"
#include "delfem2/gl_v23.h"
#include "delfem2/gl_color.h"



/////////////////////////////////////////////////////////////////////////////////////////////////////////

// display data
bool is_animation;

void draw(GLFWwindow* window)
{
  float ratio;
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  ratio = width / (float) height;
  glViewport(0, 0, width, height);
  
  ::glClearColor(1,1,1,1);
  glClear(GL_COLOR_BUFFER_BIT);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  glBegin(GL_TRIANGLES);
  glColor3f(1.f, 0.f, 0.f);
  glVertex3f(-0.6f, -0.4f, 0.f);
  glColor3f(0.f, 1.f, 0.f);
  glVertex3f(0.6f, -0.4f, 0.f);
  glColor3f(0.f, 0.f, 1.f);
  glVertex3f(0.f, 0.6f, 0.f);
  glEnd();
  
  DrawCylinder(CVector3(-1,0,0), CVector3(+1,0,0), 0.1);
  
  
  glfwSwapBuffers(window);
  glfwPollEvents();
}

static void error_callback(int error, const char* description)
{
  fputs(description, stderr);
}
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
    glfwSetWindowShouldClose(window, GL_TRUE);
  }
}

static void window_size_callback(GLFWwindow* window, int width, int height) {
  //  draw(window);
  glViewport(0, 0, width, height);
}

int main(void)
{
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()){
    exit(EXIT_FAILURE);
  }
  GLFWwindow* window = glfwCreateWindow(400, 300, "Simple example", NULL, NULL);
  if (!window){
    glfwTerminate();
    exit(EXIT_FAILURE);
  }
  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);
  glfwSetWindowSizeCallback(window, window_size_callback);
  while (!glfwWindowShouldClose(window))
  {
    draw(window);
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

