#version 120

varying vec3 LightIntensity; 
void main() 
{ 
  gl_FragColor = vec4(LightIntensity, 1.0); 
}
