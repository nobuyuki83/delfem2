#version 120

varying vec3 vLightDir; 
varying vec3 vNormal;
void main() 
{
  vLightDir = normalize( vec3( gl_LightSource[0].position ));
  vNormal = gl_Normal;
  gl_Position = ftransform();
}