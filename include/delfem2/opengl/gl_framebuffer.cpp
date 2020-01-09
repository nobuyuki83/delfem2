/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdio>

#include "glad/glad.h" // gl3.0+
#if defined(__APPLE__) && defined(__MACH__)
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif

#include "delfem2/opengl/gl_framebuffer.h"

void CFrameBufferManager::DeleteFrameBuffer(){
  if( id_framebuffer > 0 ){
    glDeleteFramebuffers(1, &id_framebuffer); // gl3.0+
    id_framebuffer = 0;
  }
  if( id_depth_render_buffer > 0  ){
    glDeleteRenderbuffers(1, &id_depth_render_buffer);
    id_depth_render_buffer = 0;
  }
  if( id_color_render_buffer > 0  ){
    glDeleteRenderbuffers(1, &id_color_render_buffer);
    id_color_render_buffer = 0;
  }
}

void CFrameBufferManager::Init
 (int width, int height,
  std::string sFormatPixelColor,
  bool isDepth)
{
  // glewInit() should be called beforehand
  this->sFormatPixelColor = sFormatPixelColor;
  DeleteFrameBuffer();
  glGenFramebuffers(1, &id_framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
  //
  glReadBuffer(GL_NONE);
  // depth
  if( isDepth ){
    glGenRenderbuffers(1, &id_depth_render_buffer);
    glBindRenderbuffer(GL_RENDERBUFFER, id_depth_render_buffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, width, height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, id_depth_render_buffer);
  }
  // color
  if( sFormatPixelColor=="4byte" || sFormatPixelColor=="4float" ){
    glGenRenderbuffers(1, &id_color_render_buffer);
    glBindRenderbuffer(GL_RENDERBUFFER, id_color_render_buffer);
    if(      sFormatPixelColor == "4byte" ){
      glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);
    }
    else if( sFormatPixelColor == "4float" ){
      glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA32F, width, height);
    }
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, id_color_render_buffer);
  }
  else{
    glDrawBuffer(GL_NONE); // not sure why do I need this...
  }
  ////
  GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER) ;
  if(status != GL_FRAMEBUFFER_COMPLETE){
    std::cout << "id_color_render_buffer: " << id_color_render_buffer << std::endl;
    std::cout << "error!: " << status << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT << std::endl;
    std::cout << GL_FRAMEBUFFER_UNSUPPORTED << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER << std::endl;
    std::cout << GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER << std::endl;
    return;
  }
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);
}

void CFrameBufferManager::Start() const{
  glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, id_depth_render_buffer);
  glBindRenderbuffer(GL_RENDERBUFFER, id_color_render_buffer);
}

void CFrameBufferManager::End() const {
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

// -------------------------------------------------------------------


