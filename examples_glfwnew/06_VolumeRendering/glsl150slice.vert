#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 mt; // texture matrix 
uniform mat4 mw; // modelview transformation
uniform mat4 mp; // projection transformation
uniform float spacing;  // space beteween slice

layout (location = 0) in vec2 pv;

out vec4 p; // mw * texture coordinate 
out vec3 t; // texture coordinate

void main()
{
  // gl_InstanceID is the index of glDrawArraysInstanced.
  p = vec4(pv, (float(gl_InstanceID) + 0.5) * spacing - 0.5, 1.0);

  // rotate and expand texture 
  t = (mat3(mt) * p.xyz) * 1.732 + 0.5;

  p = mw * p;

  gl_Position = mp * p;
}
