/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @details the file depends on OpenGL 2.x and 4.x
 * compile with -USE_GLEW if you want to use "glew", otherwise it use glad
 */

#ifndef GL24_FUNCS_H
#define GL24_FUNCS_H

#include <string>
#include <vector>

int GL24_CompileShader(const char *vert, const char* frag);

int compileShader(const std::string& str_glsl_vert,
                  int shaderType);

// compile vertex and fragment shader
// return shader program
int setUpGLSL(const std::string& str_glsl_vert,
              const std::string& str_glsl_frag);



#endif /* utility_gl_h */
