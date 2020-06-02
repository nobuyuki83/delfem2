#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 mt;
uniform mat4 mw;
uniform mat4 mp;
uniform float spacing;

layout (location = 0) in vec2 pv;

out vec4 p;
out vec3 t;

void main()
{
  p = vec4(pv, (float(gl_InstanceID) + 0.5) * spacing - 0.5, 1.0);
  t = (mat3(mt) * p.xyz) * 1.732 + 0.5;
  p = mw * p;
  gl_Position = mp * p;
}
