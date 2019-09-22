/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef GL4_FUNCS_H
#define GL4_FUNCS_H

#include <vector>
#include <assert.h>
#include <iostream> // this must be delated in future

//////////////////////////////////////////////

/**
 * @details the OpenGL ES 2.0 only accept float array. So there is no "double" version of this file
 */
int GL4_VAO_MeshTri3D
(const float* aP, int nP, int nDim,
 const unsigned int* aTri, int nTri);

class CGL4_VAO_Mesh
{
public:
  void Draw() const {
    glBindVertexArray(VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
    //glDrawArrays(GL_TRIANGLES, 0, 6);
    glDrawElements(GL_TRIANGLES, 3*nTri, GL_UNSIGNED_INT, 0);
    // glBindVertexArray(0); // no need to unbind it every time
  }
public:
  int VAO;
  int nTri;
};


const std::string glsl33vert_simplest =
"#version 330 core\n"
"layout (location = 0) in vec3 aPos;\n"
"void main()\n"
"{\n"
"   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"}\0";

const std::string glsl33vert_projection =
"#version 330 core\n"
"uniform mat4 projectionMatrix;\n"
"layout (location = 0) in vec3 aPos;\n"
"void main()\n"
"{\n"
"   gl_Position = projectionMatrix * vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"}\0";

const std::string glsl33frag =
"#version 330 core\n"
"out vec4 FragColor;\n"
"void main()\n"
"{\n"
"   FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);\n"
"}\n\0";

////////////////////////////////////////////////////

const std::string glsles3vert_simplest =
"#version 300 es\n"
"in vec3 aPos;\n"
"void main()\n"
"{\n"
"   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"}\0";

const std::string glsles3vert_projection =
"#version 300 es\n"
"uniform mat4 projectionMatrix;\n"
"in vec3 aPos;\n"
"void main()\n"
"{\n"
"   gl_Position = projectionMatrix * vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"}\0";

const std::string glsles3frag =
"#version 300 es\n"
"precision highp float;\n"
"out vec4 FragColor;\n"
"void main()\n"
"{\n"
"   FragColor = vec4(1.0, 0.5, 0.2, 1.0);\n"
"}\n\0";


#endif /* utility_glew_h */

