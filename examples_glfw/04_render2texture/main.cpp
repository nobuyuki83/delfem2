#include <iostream>
#include <cmath>
#include "delfem2/imgio.h"
#include "delfem2/noise.h"
#include "delfem2/primitive.h"

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
CShader_TriMesh shdr0;
CShader_TriMesh_Tex shdr;
delfem2::opengl::CViewer_GLFW viewer;
unsigned int targetTexture = 0;

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
  glBindTexture(GL_TEXTURE_2D , targetTexture);
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

  {
    std::vector<double> aXYZ;
    std::vector<unsigned int> aTri;
    dfm2::MeshTri3D_Torus(aXYZ, aTri, 0.8, 0.1);
    shdr0.Compile();
    shdr0.Initialize(aXYZ, aTri);
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

  const unsigned int targetTextureWidth = 256;
  const unsigned int targetTextureHeight = 256;

  {
    ::glEnable(GL_TEXTURE_2D);
    ::glActiveTexture(GL_TEXTURE0);
    // create to render to
    ::glGenTextures(1, &targetTexture);
    ::glBindTexture(GL_TEXTURE_2D, targetTexture);
    // define size and format of level 0
    ::glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                   targetTextureWidth, targetTextureHeight, 0,
                   GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    // set the filtering so we don't need mips
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }

  {
    const int nSize = 256;
    std::vector<double> aV;
    dfm2::ComputePerlin(aV,
                        nSize, nSize,
                        4, 4, 0.8);
    assert( aV.size() == nSize*nSize );
    std::vector<unsigned char> image(nSize*nSize*3);
    for(unsigned int i=0;i<aV.size();i++){
      int ival = int( aV[i]*50 + 128 );
      if( ival < 0 ){ ival = 0; }
      if( ival >= 256 ){ ival = 255; }
      image[i*3+0] = ival;
      image[i*3+1] = ival;
      image[i*3+2] = ival;
    }
    // -------
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D , targetTexture);
    glTexImage2D(GL_TEXTURE_2D , 0 , GL_RGB , targetTextureWidth, targetTextureHeight,
                 0 , GL_RGB , GL_UNSIGNED_BYTE , image.data() );
    glGenerateMipmap(GL_TEXTURE_2D);
  }

  {
    // Create and bind the framebuffer
    unsigned int fb;
    ::glGenFramebuffers(1,&fb);
    ::glBindFramebuffer(GL_FRAMEBUFFER, fb);

    // attach the texture as the first color attachment
    ::glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
        targetTexture, 0);

    // Always check that our framebuffer is ok
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      return 0;

    int viewport[4];
    ::glGetIntegerv(GL_VIEWPORT,viewport);

    ::glViewport(0,0,targetTextureWidth,targetTextureHeight);
//    ::glClearColor(0.8, 0.8, 0.8, 1.0);
//    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    float mP[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    float mMV[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    shdr0.Draw(mP,mMV);

    ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
    ::glViewport(0,0,viewport[2],viewport[3]);
  }

  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

