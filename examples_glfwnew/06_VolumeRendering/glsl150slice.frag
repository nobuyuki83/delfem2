#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

const vec4 lamb = vec4(0.4, 0.4, 0.4, 1.0);
const vec4 ldiff = vec4(0.8, 0.8, 0.8, 1.0);
const vec4 lspec = vec4(0.8, 0.8, 0.8, 1.0);
const vec4 pl = vec4(0.0, 0.5, 1.0, 0.0);
const vec4 kamb = vec4(0.8, 0.8, 0.8, 1.0);
const vec4 kdiff = vec4(0.8, 0.6, 0.6, 1.0);
const vec4 kspec = vec4(0.4, 0.4, 0.4, 1.0);
const float kshi = 50.0;

uniform sampler3D volume;
uniform mat4 mt;
uniform float threshold;

in vec4 p;
in vec3 t;

layout (location = 0) out vec4 fc;

void main()
{
  float v = texture(volume, t).r - threshold;
  if (v <= 0.0) discard;

  vec3 g = vec3(
    textureOffset(volume, t, ivec3(-1, 0, 0)).r - textureOffset(volume, t, ivec3(1, 0, 0)).r,
    textureOffset(volume, t, ivec3(0, -1, 0)).r - textureOffset(volume, t, ivec3(0, 1, 0)).r,
    textureOffset(volume, t, ivec3(0, 0, -1)).r - textureOffset(volume, t, ivec3(0, 0, 1)).r
  );

  vec3 l = normalize((pl * p.w - p * pl.w).xyz);
  vec3 n = normalize(g * mat3(mt));
  vec3 h = normalize(l - normalize(p.xyz));
  vec4 idiff = max(dot(n, l), 0.0) * kdiff * ldiff + kamb * lamb;
  vec4 ispec = pow(max(dot(n, h), 0.0), kshi) * kspec * lspec;
  fc = vec4((idiff + ispec).rgb, v);
//  fc = vec4(normalize(g) * 0.5 + 0.5, v); // draw normal

}
