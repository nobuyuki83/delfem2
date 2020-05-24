#version 120

varying vec2 vTexCoord;

uniform sampler2D TexOrg;
uniform sampler2D TexNrm;
uniform int nTexWidth;
uniform int nTexHeight;

void main()
{
  vec3 clr = texture2D(TexOrg, vTexCoord ).rgb;

  vec3 clr00 = texture2D(TexNrm, vTexCoord + vec2(-1.0/nTexWidth, -1.0/nTexHeight) ).rgb;
  vec3 clr01 = texture2D(TexNrm, vTexCoord + vec2(-1.0/nTexWidth, +0.0/nTexHeight) ).rgb;
  vec3 clr02 = texture2D(TexNrm, vTexCoord + vec2(-1.0/nTexWidth, +1.0/nTexHeight) ).rgb;
  vec3 clr10 = texture2D(TexNrm, vTexCoord + vec2(+0.0/nTexWidth, -1.0/nTexHeight) ).rgb;
  vec3 clr11 = texture2D(TexNrm, vTexCoord + vec2(+0.0/nTexWidth, +0.0/nTexHeight) ).rgb;
  vec3 clr12 = texture2D(TexNrm, vTexCoord + vec2(+0.0/nTexWidth, +1.0/nTexHeight) ).rgb;
  vec3 clr20 = texture2D(TexNrm, vTexCoord + vec2(+1.0/nTexWidth, -1.0/nTexHeight) ).rgb;
  vec3 clr21 = texture2D(TexNrm, vTexCoord + vec2(+1.0/nTexWidth, +0.0/nTexHeight) ).rgb;
  vec3 clr22 = texture2D(TexNrm, vTexCoord + vec2(+1.0/nTexWidth, +1.0/nTexHeight) ).rgb;

  vec3 lap = 8*clr11-clr00-clr01-clr02-clr10-clr12-clr20-clr21-clr22;
  vec3 tmpV = vec3(1.0);
  float d = abs(dot(lap,tmpV))/3.0;
  if( d > 1){ d = 1; }
////  gl_FragColor = vec4(clr11, 1.0);
////  gl_FragColor = vec4(clr, 1.0);
  gl_FragColor = vec4((1.0-d)*clr, 1.0);
//   gl_FragColor = vec4(1-d,1-d,1-d, 1.0);
//  gl_FragColor = vec4(lap, 1.0);
}
