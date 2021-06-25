#version 120

varying vec3 vNormal;

void main(){ 
  gl_Position = ftransform();
  vNormal   = vec3(gl_NormalMatrix    * gl_Normal);
}
