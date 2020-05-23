#version 120

varying vec2 vTexCoord;

void main()
{
  gl_Position = ftransform();
  vTexCoord = gl_MultiTexCoord0.xy;
}
