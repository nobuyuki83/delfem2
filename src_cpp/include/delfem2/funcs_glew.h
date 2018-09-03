#ifndef FUNCS_GLEW_H
#define FUNCS_GLEW_H

#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#else
#include <GL/glew.h>
#endif

static int compileShader
(const std::string& str_glsl_vert,
 int shaderType)
{
  int id_shader = glCreateShader(shaderType);
  const char *vfile = str_glsl_vert.c_str();
  glShaderSource(id_shader, 1, &vfile, NULL);
  glCompileShader(id_shader); // compile the code
  
  {
    GLint res;
    glGetShaderiv(id_shader, GL_COMPILE_STATUS, &res);
    if (res==GL_FALSE){
      if (shaderType==GL_VERTEX_SHADER){
        std::cout<<"compile vertex shader failed"<<std::endl;
      }
      else if(shaderType==GL_FRAGMENT_SHADER){
        std::cout<<"compile fragment shader failed"<<std::endl;
      }
    }
  }
  return id_shader;
}

// compile vertex and fragment shader
// return shader program
static int setUpGLSL
(const std::string& str_glsl_vert,
 const std::string& str_glsl_frag)
{
  int vShaderId = compileShader(str_glsl_vert, GL_VERTEX_SHADER);
  int fShaderId = compileShader(str_glsl_frag, GL_FRAGMENT_SHADER);
  
  int id_program = glCreateProgram();
  glAttachShader(id_program,vShaderId);
  glAttachShader(id_program,fShaderId);
  
  GLint linked;
  glLinkProgram(id_program);
  glGetProgramiv(id_program, GL_LINK_STATUS, &linked);
  if(linked == GL_FALSE)
  {
    std::cerr << "Link Err.\n";
  }
  return id_program;
}

class CFrameBufferManager
{
public:
  CFrameBufferManager(){
    id_framebuffer = -1;
    id_depth_texture = -1;
  }
  CFrameBufferManager(const std::vector<int>& winSize){
    this->Init(winSize[0],winSize[1]);
  }
  void Init(int width, int height){
    this->width = width;
    this->height = height;
    this->SetFrameBufferSize(width,height);
  }
  void DeleteFrameBuffer(){
    if( id_framebuffer > 0 ){
      glDeleteFramebuffers(1, &id_framebuffer);
      id_framebuffer = 0;
    }
    // TODO: delete depth_texture here
    if( id_depth_render_buffer > 0  ){
      glDeleteRenderbuffersEXT(1, &id_depth_render_buffer);
      id_depth_render_buffer = 0;
    }
  }
  void SetFrameBufferSize(int width, int height){
    DeleteFrameBuffer();
    glGenFramebuffers(1, &id_framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
    ////
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);
    ////
    glGenRenderbuffers(1, &id_depth_render_buffer);
    glBindRenderbuffer(GL_RENDERBUFFER, id_depth_render_buffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, width, height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, id_depth_render_buffer);
    ////
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER) ;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    if(status != GL_FRAMEBUFFER_COMPLETE){
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
  }
  void Start() const{
    glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, id_depth_render_buffer);
  }
  void End() const {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }
public:
  unsigned int id_framebuffer;
  unsigned int id_depth_texture;
  unsigned int id_depth_render_buffer;
  int width;
  int height;
};

#endif /* utility_glew_h */

