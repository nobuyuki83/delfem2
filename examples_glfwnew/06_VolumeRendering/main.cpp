
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include "delfem2/noise.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "delfem2/opengl/gl_funcs.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;
 
GLuint MakeVAO_Slice()
{
  GLuint vao; glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  //
  GLuint vbo; glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  {
    const GLfloat p[8] = {
      -0.5f, -0.5f,
      +0.5f, -0.5f,
      +0.5f, +0.5f,
      -0.5f, +0.5f };
    glBufferData(GL_ARRAY_BUFFER, sizeof p, p, GL_STATIC_DRAW);
  }
   // define an array of generic vertex attrib data
  glVertexAttribPointer(0, // index
                        2, // size
                        GL_FLOAT, // type
                        GL_FALSE, // normalized
                        0, // stride
                        0); // pointer
  glEnableVertexAttribArray(0); // enable the previous vertex attrib array
  return vao;
}



// making 3D volumetric texture
GLuint MakeTex_Volume(GLint nW, GLint nH, GLint nD)
{
  std::vector<int> aP;
  {
    aP.resize(256);
    for(int i=0;i<256;++i){ aP[i]=i; }
    delfem2::Shuffle(aP);
    aP.resize(512);
    for(int i=0;i<256;++i){ aP[256+i]=i; }
  }
  
  std::vector<double> aGrad;
  aGrad.push_back(-1); aGrad.push_back(-1); aGrad.push_back(+0);
  aGrad.push_back(-1); aGrad.push_back(+1); aGrad.push_back(+0);
  aGrad.push_back(+1); aGrad.push_back(-1); aGrad.push_back(+0);
  aGrad.push_back(+1); aGrad.push_back(+1); aGrad.push_back(+0);
  aGrad.push_back(+0); aGrad.push_back(-1); aGrad.push_back(-1);
  aGrad.push_back(+0); aGrad.push_back(-1); aGrad.push_back(+1);
  aGrad.push_back(+0); aGrad.push_back(+1); aGrad.push_back(-1);
  aGrad.push_back(+0); aGrad.push_back(+1); aGrad.push_back(+1);
  aGrad.push_back(-1); aGrad.push_back(+0); aGrad.push_back(-1);
  aGrad.push_back(-1); aGrad.push_back(+0); aGrad.push_back(+1);
  aGrad.push_back(+1); aGrad.push_back(+0); aGrad.push_back(-1);
  aGrad.push_back(+1); aGrad.push_back(+0); aGrad.push_back(+1);
  
  std::vector<GLubyte> aV; // temporary buffer
  aV.resize(nH*nW*nD*4);
  int nrep = 4;
  for(int id=0;id<nD;++id){
    for(int ih=0;ih<nH;++ih){
      for(int iw=0;iw<nW;++iw){
        double x = (double)iw/nH*nrep;
        double y = (double)ih/nW*nrep;
        double z = (double)id/nD*nrep;
        double v = delfem2::noise_perlin_3d_oct(x,y,z,nrep, 4,0.5, aGrad,aP);
        //        double v = noise_perlin_3d(x,y,z, aGrad,aP);
        double v0 = v*128+128;
        if( v0 < 0   ){ v0 =   0; }
        if( v0 > 255 ){ v0 = 255; }
        auto ucv = (unsigned char)v0;
        aV[id*nW*nH+ih*nW+iw] = ucv;
      }
    }
  }
  
  // generate 3D texture and bind it
  GLuint tex; glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_3D, tex);
  
  // send 3D texture data to GPU
  glTexImage3D(GL_TEXTURE_3D,
               0, GL_R8, nW, nH, nD, 0,
               GL_RED, GL_UNSIGNED_BYTE, aV.data());

  // linear interplation
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  // use border color outside texture range
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);

  { // specify border color
    const GLfloat border_color[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    glTexParameterfv(GL_TEXTURE_3D, GL_TEXTURE_BORDER_COLOR, border_color);
  }

  return tex;
}

std::string LoadFile
 (const std::string& fname)
{
  std::ifstream inputFile1(fname.c_str());
  std::istreambuf_iterator<char> vdataBegin(inputFile1);
  std::istreambuf_iterator<char> vdataEnd;
  return std::string(vdataBegin,vdataEnd);
}

GLuint MakeTex_Gradient(GLuint idVAO,
                        GLuint idTexVolume,
                        GLint texWidth,
                        GLint texHeight,
                        GLint texDepth)
{
  GLuint idTex; glGenTextures(1, &idTex);
  glBindTexture(GL_TEXTURE_3D, idTex);
  
  glTexImage3D(GL_TEXTURE_3D,
               0, GL_RGBA, texWidth, texHeight, texDepth, 0,
               GL_RGBA, GL_UNSIGNED_BYTE, 0);
  
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  
  glBindTexture(GL_TEXTURE_3D, 0);
  
  // get viewport
  GLint viewport[4]; glGetIntegerv(GL_VIEWPORT, viewport);
  
  // set up viewport
  glViewport(0, 0, texWidth, texHeight);
  
  GLuint idFrameBuffer; glGenFramebuffers(1, &idFrameBuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, idFrameBuffer);
  
  static const GLenum bufs[] = { GL_COLOR_ATTACHMENT0 };
  glDrawBuffers(1, bufs);
  
  glDepthMask(GL_FALSE); // don't use depth buffer
  glDisable(GL_BLEND); // disable alpha blending
  
  GLuint idProgramGradient = 0;
  {
    std::string sgradv = LoadFile(std::string(PATH_SOURCE_DIR)+std::string("/glsl150gradient.vert"));
    std::string sgradf = LoadFile(std::string(PATH_SOURCE_DIR)+std::string("/glsl150gradient.frag"));
    idProgramGradient = dfm2::opengl::GL24_CompileShader(sgradv.c_str(),sgradf.c_str());
  }
  GLint locTex = glGetUniformLocation(idProgramGradient, "tex");
  GLint locZ = glGetUniformLocation(idProgramGradient, "z");
  
  glUseProgram(idProgramGradient);
  glUniform1i(locTex, 0);
  
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_3D, idTexVolume);
  glBindVertexArray(idVAO);
  for (GLint z = 0; z < texDepth; ++z)
  {
    glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, idTex, 0, z);
    glUniform1f(locZ, (static_cast<GLfloat>(z) + 0.5f) / static_cast<GLfloat>(texDepth));
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
  }
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  
  glDrawBuffer(GL_BACK);
  glDepthMask(GL_TRUE);
  
  // restore the original viewport
  glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
  
  return idTex;
}

int main(int argc, const char * argv[])
{
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_newGL();
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  // making VAO for slice
  const GLuint idVaoSlice = MakeVAO_Slice();

  const GLsizei texWidth = 32;
  const GLsizei texHeight = 32;
  const GLsizei texDepth = 32;
  
  // making 3D texture
  const GLuint idTexVolume = MakeTex_Volume(texWidth, texHeight, texDepth);

  // making grdient texture
  const GLuint idTexGradient = MakeTex_Gradient(idVaoSlice, idTexVolume, texWidth, texHeight, texDepth);

  GLuint idProgramSlice = 0;
  { // compile slicing shader
    const std::string sslicev = LoadFile(std::string(PATH_SOURCE_DIR)+"/glsl150slice.vert");
    const std::string sslicef = LoadFile(std::string(PATH_SOURCE_DIR)+"/glsl150slice.frag");
    idProgramSlice = dfm2::opengl::GL24_CompileShader(sslicev.c_str(),sslicef.c_str());
  }
  
  const GLint locGradient = glGetUniformLocation(idProgramSlice, "gradient");
  const GLint locVolume = glGetUniformLocation(idProgramSlice, "volume");
  const GLint locMt = glGetUniformLocation(idProgramSlice, "mt");
  const GLint locMw = glGetUniformLocation(idProgramSlice, "mw");
  const GLint locMp = glGetUniformLocation(idProgramSlice, "mp");
  const GLint locSpacing = glGetUniformLocation(idProgramSlice, "spacing");
  const GLint locThreshold = glGetUniformLocation(idProgramSlice, "threshold");
  
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);

  // set up alpha blending mode
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  { // background color
    const GLfloat bgcolor[4] = { 0.2f, 0.3f, 0.4f, 0.0f };
    glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], bgcolor[3]);
  }

  while (true)
  { // draw loop
    ::glClearColor(0.8, 1.0, 1.0, 1.0);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_POLYGON_OFFSET_FILL );
    ::glPolygonOffset( 1.1f, 4.0f );
    
    float mP[16], mMV[16];
    viewer.nav.Matrix_MVP(mMV, mP, viewer.window);
    
    const unsigned int nslice = 256;

    glUseProgram(idProgramSlice);
    glUniformMatrix4fv(locMt, 1, GL_TRUE, mMV);
    glUniform1f(locSpacing, 1.0f / static_cast<GLfloat>(nslice - 1));
    glUniform1f(locThreshold, 0.5);

    glUniform1i(locVolume, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, idTexVolume);
    glUniform1i(locGradient, 1);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, idTexGradient);

    glBindVertexArray(idVaoSlice); // bind VAO
    glUniformMatrix4fv(locMw, 1, GL_FALSE, mMV);
    glUniformMatrix4fv(locMp, 1, GL_FALSE, mP);
    glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, nslice); // draw
    
    viewer.DrawEnd_oldGL();
    
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
  return 0;
}
