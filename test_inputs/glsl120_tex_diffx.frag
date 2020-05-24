#version 120

varying vec2 vTexCoord;

uniform sampler2D ourTexture;
uniform int nTexWidth;
uniform int nTexHeight;

void main()
{
  vec3 clr = texture2D(ourTexture, vTexCoord ).rgb;

  vec2 tc1 = vTexCoord + vec2(1.0/nTexWidth, 0.0/nTexHeight);
  vec3 clr1 = texture2D(ourTexture, tc1 ).rgb;

  vec3 dclr = abs(clr-clr1);
  gl_FragColor = vec4(dclr, 1.0);
}
