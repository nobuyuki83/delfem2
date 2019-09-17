#include<stdio.h>
#include<stdlib.h>
#ifdef EMSCRIPTEN
#include<emscripten/emscripten.h>
#define GLFW_INCLUDE_ES3
#endif
#include <GLFW/glfw3.h>
#include <iostream>

int windowWidth;
int windowHeight;

void window_size_callback(GLFWwindow* window, int width, int height)
{
    printf("window_size_callback received width: %i, height: %i\n", width, height);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int modifier)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE)
        glfwSetWindowShouldClose(window, 1);

    if (key == GLFW_KEY_ENTER)
        std::cout << "Enter was hit\n" << std::endl;
}

static void wmbutcb(GLFWwindow *window, int button, int action, int mods) {
    //assert(window != NULL); (void)button; (void)action; (void)mods;
    printf("Mouse buttion! \n");
}

GLFWwindow *window;

void do_frame(){
    static int a = 0;
    printf("Fc: %d \n", ++a);
    glClearColor(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glfwSwapBuffers(window);
    glfwPollEvents();
}

int main(int argc, char **argv) {


if (glfwInit()!=GL_TRUE) {
    printf("glfwInit() failed\n");
    glfwTerminate();
} else {
    printf("glfwInit() success\n");
    window = glfwCreateWindow(1280, 512, "GLFW test", NULL, NULL);
    if (!window){
        printf("glfwCreateWindow() failed\n");
        glfwTerminate();
    } else {
      printf("glfwCreateWindow() success\n");
      glfwMakeContextCurrent(window);
      glfwGetFramebufferSize(window, &windowWidth, &windowHeight);
      glfwSetWindowSizeCallback(window, window_size_callback);
      glfwSetMouseButtonCallback(window, wmbutcb);
      glfwSetKeyCallback(window, key_callback);
      glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
#ifdef EMSCRIPTEN
      emscripten_set_main_loop(do_frame, 0, 1);
#else
      while (!glfwWindowShouldClose(window))
      {
          do_frame();
      }
#endif
    }
}

    glfwTerminate();
return EXIT_SUCCESS;
}
