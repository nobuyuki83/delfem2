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

#include "delfem2/cad2d.h"

#include "delfem2/gl24_funcs.h"
#include "delfem2/gl4_funcs.h"
#include "delfem2/gl4_v23dtricad.h"
#include "../glfw_funcs.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////

CNav3D_GLFW nav;
CCad2D cad;
CShader_CCad2D shdr_cad;

void draw(GLFWwindow* window)
{
  float asp;
  {
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    asp = width / (float) height;
//    std::cout << width << " " << height << " " << asp << std::endl;
  }
  
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );

  float mP[16]; nav.camera.Affine4f_Projection(mP, asp, 10);
  float mMV[16]; nav.camera.Affine4f_ModelView(mMV);
  shdr_cad.Draw(mP,mMV);
  
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

void callback_mouse_button(GLFWwindow* window, int button, int action, int mods)
{
  nav.Mouse(window,button,action,mods);

}

void callback_cursor_position(GLFWwindow* window, double xpos, double ypos)
{
  nav.Motion(window,xpos,ypos);
}

void callback_scroll(GLFWwindow* window, double xoffset, double yoffset)
{
  nav.camera.scale *= pow(1.01,yoffset);
}


int main(void)
{
  GLFWwindow* window = myGLFW_OpenWindow(800,600);
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, callback_resize);
  glfwSetKeyCallback(            window, callback_key);
  glfwSetMouseButtonCallback(    window, callback_mouse_button);
  glfwSetCursorPosCallback(      window, callback_cursor_position);
  glfwSetScrollCallback(         window, callback_scroll);
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  {
    {
      std::vector<double> aXY = {-1,-1, +1,-1, +1,+1, -1,+1};
      cad.AddPolygon(aXY);
    }
    shdr_cad.MakeBuffer(cad);
    shdr_cad.Compile();
    {
      std::vector<int> aFlgPnt, aFlgTri;
      CMeshDynTri2D dmsh;
      CMesher_Cad2D mesher;
      mesher.edge_length = 0.05;
      mesher.Meshing(dmsh, cad);
    }
  }
  
  nav.camera.view_height = 1.5;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  while (!glfwWindowShouldClose(window)) { draw(window); }
#endif
  
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

