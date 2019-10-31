#include <iostream>
#include <math.h>

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/gl4_funcs.h"
#include "delfem2/gl24_funcs.h"
#include "delfem2/glfw_cam.h"

// -------------------------------------------

int shaderProgram;
unsigned int VAO;
unsigned int VBO_Tri;
unsigned int EBO_Tri;

void draw(GLFWwindow* window)
{
  float ratio;
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  ratio = width / (float) height;
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



int main(void)
{
  GLFWwindow* window = myGLFW_OpenWindow(800,600);
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, callback_resize);
  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, callback_key);
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

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
    GL4_VAO_Pos(VAO, VBO_Tri,
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
    shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
                                        glslvrt_simplest).c_str(),
                                       (std::string("#version 300 es\n")+
                                        std::string("precision highp float;\n")+
                                        glslfrg_simplest).c_str());
#else
    shaderProgram = GL24_CompileShader(("#version 330 core\n"+glslvrt_simplest).c_str(),
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

