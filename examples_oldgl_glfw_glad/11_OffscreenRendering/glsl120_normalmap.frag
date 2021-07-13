#version 120

varying vec3 normal;
void main()
{
  gl_FragColor = vec4(0.5*normal+0.5,1);
}
