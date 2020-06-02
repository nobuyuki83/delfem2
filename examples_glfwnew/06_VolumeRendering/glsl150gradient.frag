#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform sampler3D tex;

in vec3 t;

layout (location = 0) out vec4 fc;

void main(void)
{
  vec3 g = vec3(
    textureOffset(tex, t, ivec3(-1, 0, 0)).r - textureOffset(tex, t, ivec3(1, 0, 0)).r,
    textureOffset(tex, t, ivec3(0, -1, 0)).r - textureOffset(tex, t, ivec3(0, 1, 0)).r,
    textureOffset(tex, t, ivec3(0, 0, -1)).r - textureOffset(tex, t, ivec3(0, 0, 1)).r
  );

  fc = vec4(normalize(g) * 0.5 + 0.5, length(g));
}
