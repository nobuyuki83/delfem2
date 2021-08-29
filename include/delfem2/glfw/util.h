#include <GLFW/glfw3.h>
#include <cstdio>
#include <cstdlib>

namespace delfem2{
namespace glfw{

static void callback_error(
    [[maybe_unused]] int error,
    const char *description) {
  fputs(description, stderr);
}

void InitGLOld()
{
  glfwSetErrorCallback(callback_error);
  if (!glfwInit()){
    exit(EXIT_FAILURE);
  }
  // shader version #120
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
}

void InitGLNew(){
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
}

}
}
