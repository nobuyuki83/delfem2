#include <iostream>
#include <cmath>
#include "delfem2/imgio.h"

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/opengl/gl24_funcs.h"
#include "delfem2/opengl/gl24_tex.h"
#include "delfem2/opengl/gl4_mshcolor.h"
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/glfw_cam.h"

namespace dfm2 = delfem2;

// ---------------------------
CShader_TriMesh_Tex shdr;
delfem2::opengl::CViewer_GLFW viewer;
GLuint m_texName = -1;
//GLuint sampler_obj= 0;


// ---------------------------


void draw(GLFWwindow* window)
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//  ::glEnable(GL_DEPTH_TEST);
//  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );

  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0); // activate the texture unit first before binding texture
  glBindTexture(GL_TEXTURE_2D , m_texName);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );

  float mP[16], mMV[16];
  viewer.nav.Matrix_MVP(mMV, mP, window);
  shdr.Draw(mP, mMV);
  
  viewer.DrawEnd_oldGL();
}

int main()
{
  viewer.Init_newGL();
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  shdr.Compile();
  {
    std::vector<double> aPos3d = {
        -0.5, -0.5, 0.0,
        +0.5, -0.5, 0.0,
        +0.5, +0.5, 0.0,
        -0.5, +0.5, 0.0
    };
    std::vector<unsigned int> aTri = {
        0,1,2,
        0,2,3,
    };
    std::vector<double> aTex2d = {
        0.0, 0.0,
        1.0, 0.0,
        1.0, 1.0,
        0.0, 1.0
    };
    shdr.Initialize(aPos3d, aTri, aTex2d);
  }

  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0);
  const std::string path = std::string(PATH_INPUT_DIR)+"/dep.ppm";
  {
    std::cout << path << std::endl;
    unsigned int w, h;
    std::vector<unsigned char> image;
    dfm2::LoadImage_PPMAscii(w, h, image,
        path);
    assert(image.size() == w * h * 3);
    m_texName = dfm2::opengl::SetTexture_RGB(w,h,image);
  }
  glGenerateMipmap(GL_TEXTURE_2D);

  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  
#ifdef EMSCRIPTEN
  std::cout << "I know I cannot open file in Emsripten." << std::endl;
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

