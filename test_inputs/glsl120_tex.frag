#version 120

varying vec2 vTexCoord;

uniform sampler2D ourTexture;
uniform int nTexWidth;
uniform int nTexHeight;

void main()
{
  vec3 clr = texture2D(ourTexture, vTexCoord ).rgb;
  gl_FragColor = vec4(clr, 1.0);
}
