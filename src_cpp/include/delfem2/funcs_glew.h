#ifndef FUNCS_GLEW_H
#define FUNCS_GLEW_H

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

#endif /* utility_glew_h */

