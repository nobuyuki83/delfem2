#version 120

varying vec3 vNormal;
varying vec3 vPosition;

void main(){ 
  gl_Position = ftransform();
  vNormal   = vec3(gl_NormalMatrix    * gl_Normal);
//  vPosition = vec3(gl_ModelViewMatrix * gl_Vertex);
//  vPosition = vec3(gl_ModelViewMatrix * vec4(0,0,1,0));
  vPosition = vec3(0,0,1);
}
