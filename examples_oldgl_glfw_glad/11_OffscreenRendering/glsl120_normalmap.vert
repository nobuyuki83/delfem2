#version 120

varying vec3 normal;
void main() 
{
  gl_Position = ftransform();
  normal = vec3(gl_NormalMatrix*gl_Normal);
}
