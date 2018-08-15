varying vec3 vNormal; 
varying vec3 vPosition;

void main(){ 
  gl_Position = ftransform();
  vNormal   = vec3(gl_NormalMatrix    * gl_Normal);
  vPosition = vec3(gl_ModelViewMatrix * gl_Vertex);
}