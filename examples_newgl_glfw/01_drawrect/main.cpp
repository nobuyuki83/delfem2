#include <iostream>
#include <cmath>

#if defined(_MSC_VER)
#  include <windows.h>
#endif

#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#else
#  include <glad/glad.h> // emscripten and glad didn't work together
#endif
#include "delfem2/opengl/funcs.h"
#include "delfem2/opengl/new/funcs.h"
#include <GLFW/glfw3.h>


namespace dfm2 = delfem2;

static void callback_error(
    [[maybe_unused]] int error,
    const char* description)
{
  fputs(description, stderr);
}

static GLFWwindow* myGLFW_OpenWindow(
    const unsigned int SCR_WIDTH,
    const unsigned int SCR_HEIGHT)
{
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
  GLFWwindow* window = glfwCreateWindow( SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
  if (window == NULL)
  {
//    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return 0;
  }
  return window;
}



// -------------------------------------------

int shaderProgram;
unsigned int VAO;
unsigned int VBO_Tri;
unsigned int EBO_Tri;

void draw(GLFWwindow* window)
{
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  glViewport(0, 0, width, height);

  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(shaderProgram);
  glBindVertexArray(VAO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri);
  glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

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
  GLFWwindow* window = myGLFW_OpenWindow(800,600);
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, callback_resize);
  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, callback_key);

  // glad: load all OpenGL function pointers
  // ---------------------------------------
#ifndef EMSCRIPTEN
  // emscripten and glad didn't work well together
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

  {
    float vertices[] = {
        0.5f,  0.5f, 0.0f,  // top right
        0.5f, -0.5f, 0.0f,  // bottom right
        -0.5f, -0.5f, 0.0f,  // bottom left
        -0.5f,  0.5f, 0.0f   // top left
    };
    unsigned int indices[] = {  // note that we start from 0!
        0, 1, 3,  // first Triangle
        1, 2, 3   // second Triangle
    };
    dfm2::opengl::GL4_VAO_Pos(VAO, VBO_Tri,
                              vertices,4,3);
    {
      glBindVertexArray(VAO); // opengl4
      glGenBuffers(1, &EBO_Tri);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*6, indices, GL_STATIC_DRAW);
    }
  }

  {
    const std::string glslvrt_simplest =
        "layout (location = 0) in vec3 aPos;\n"
        "void main()\n"
        "{\n"
        "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
        "}\0";

    const std::string glslfrg_simplest =
        "out vec4 FragColor;\n"
        "void main()\n"
        "{\n"
        "  FragColor = vec4(1.0f, 0.5f, 0.5f, 1.0f);\n"
        "}\n\0";
#ifdef EMSCRIPTEN
    shaderProgram = dfm2::opengl::GL24_CompileShader((std::string("#version 300 es\n")+
                                                      glslvrt_simplest).c_str(),
                                                     (std::string("#version 300 es\n")+
                                                      std::string("precision highp float;\n")+
                                                      glslfrg_simplest).c_str());
#else
    shaderProgram = dfm2::opengl::GL24_CompileShader(("#version 330 core\n"+glslvrt_simplest).c_str(),
                                                     ("#version 330 core\n"+glslfrg_simplest).c_str());
#endif
  }

#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  while (!glfwWindowShouldClose(window)) { draw(window); }
#endif

  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

