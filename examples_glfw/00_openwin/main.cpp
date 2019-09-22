/**
 * @brief just open a window (should be compiled with both g++ and emscripten)
 * @details this file need to be independent as much as possible
 */

#include<stdio.h>
#include<stdlib.h>
#ifdef EMSCRIPTEN
#include <emscripten/emscripten.h>
#define GLFW_INCLUDE_ES3
#endif
#include <GLFW/glfw3.h>
#include <iostream>

void callback_resize(GLFWwindow* window, int width, int height)
{
  printf("window_size_callback received width: %i, height: %i\n", width, height);
}

void callback_key(GLFWwindow* window, int key, int scancode, int action, int modifier)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE)
    glfwSetWindowShouldClose(window, 1);
  
  if (key == GLFW_KEY_ENTER)
    std::cout << "Enter was hit\n" << std::endl;
}

void callback_mouse(GLFWwindow *window, int button, int action, int mods) {
  //assert(window != NULL); (void)button; (void)action; (void)mods;
  printf("Mouse buttion! \n");
}

void do_frame(GLFWwindow *window)
{
  /*
    static int a = 0;
    printf("Fc: %d \n", ++a);
    glClearColor(rand() / (float)RAND_MAX,
                 rand() / (float)RAND_MAX,
                 rand() / (float)RAND_MAX,
                 1.0f);
   */
  glClearColor(0.8f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  glfwSwapBuffers(window);
  glfwPollEvents();
}

int main(int argc, char **argv) {
  if (glfwInit()!=GL_TRUE) {
    printf("glfwInit() failed\n");
    glfwTerminate();
    return 0;
  }
  printf("glfwInit() success\n");
  GLFWwindow* window = glfwCreateWindow(512, 512, "GLFW test", NULL, NULL);
  
  if (!window){
    printf("glfwCreateWindow() failed\n");
    glfwTerminate();
    return 0;
  }
  
  printf("glfwCreateWindow() success\n");
  glfwMakeContextCurrent(window);
  int windowWidth;
  int windowHeight;
  glfwGetFramebufferSize(window, &windowWidth, &windowHeight);
  glfwSetWindowSizeCallback(window, callback_resize);
  glfwSetMouseButtonCallback(window, callback_mouse);
  glfwSetKeyCallback(window, callback_key);
  glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
#ifdef EMSCRIPTEN
//  emscripten_set_main_loop(do_frame, 0, 1);
  emscripten_set_main_loop_arg((em_arg_callback_func) do_frame, window, 60, 1);
#else
  while (!glfwWindowShouldClose(window)){ do_frame(window); }
#endif
  
  glfwTerminate();
  return EXIT_SUCCESS;
}
