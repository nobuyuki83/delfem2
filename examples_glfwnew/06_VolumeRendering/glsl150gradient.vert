#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform float z;

layout (location = 0) in vec2 pv;

out vec3 t;

void main(void)
{
  t = vec3(pv + 0.5, z);
  gl_Position = vec4(t * 2.0 - 1.0, 1.0);
}
