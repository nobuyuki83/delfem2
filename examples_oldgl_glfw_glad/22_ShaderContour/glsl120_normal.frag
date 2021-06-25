#version 120

varying vec3 vNormal;

void main() 
{
  vec3 norm = normalize(vNormal);
  vec3 rgb = norm*0.5+0.5;
	gl_FragColor = vec4(rgb,1.0);
}
