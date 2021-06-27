#version 330 core

uniform vec2 resolution;
uniform mat4 mMVPinv;
uniform mat4 mMV;
out vec4 FragColor;

float squareRayTest(
  vec3 org,
  vec3 dir,
  vec3 center,
  float hlen,
  vec3 n0,
  vec3 t1,
  vec3 t2 )
{
  float a1 = dot(org-center,n0)*dot(dir,t1);
  float b1 = dot(org-center,t1)*dot(dir,n0);
  float a2 = dot(org-center,n0)*dot(dir,t2);
  float b2 = dot(org-center,t2)*dot(dir,n0);
  float c = hlen*abs(dot(dir,n0));
  bool is1 = abs(a1-b1) < c;
  bool is2 = abs(a2-b2) < c;
  if( !is1 || !is2 ){ return -1.0; }
  return dot(center-org,n0)/dot(dir,n0);
}

void cubeRayTest(
  out float tmin,
  out float tmax,
  vec3 org,
  vec3 dir,
  float hlen)
{
  tmin = -1;
  tmax = -1;
  {
    float tz0 = squareRayTest(org, dir, vec3(0, 0, -hlen), hlen, vec3(0, 0, -1), vec3(1, 0, 0), vec3(0, 1, 0));
    if (tz0 > 0 && ( tmin<0 || tz0 < tmin ) ){ tmin = tz0; }
    tmax = max(tmax, tz0);
  }
  {
    float tz1 = squareRayTest( org, dir, vec3(0,0,+hlen), hlen, vec3(0,0,+1), vec3(1,0,0), vec3(0,1,0));
    if (tz1 > 0 && ( tmin<0 || tz1 < tmin ) ){ tmin = tz1; }
    tmax = max(tmax,tz1);
  }
  {
    float tx0 = squareRayTest(org, dir, vec3(-hlen, 0, 0), hlen, vec3(-1, 0, 0), vec3(0, 0, 1), vec3(0, 1, 0));
    if (tx0 > 0 && ( tmin<0 || tx0 < tmin ) ){ tmin = tx0; }
    tmax = max(tmax, tx0);
  }
  {
    float tx1 = squareRayTest(org, dir, vec3(+hlen, 0, 0), hlen, vec3(+1, 0, 0), vec3(0, 0, 1), vec3(0, 1, 0));
    if (tx1 > 0 && ( tmin<0 || tx1 < tmin ) ){ tmin = tx1; }
    tmax = max(tmax, tx1);
  }
  {
    float ty0 = squareRayTest(org, dir, vec3(0, -hlen, 0), hlen, vec3(0, -1, 0), vec3(0, 0, 1), vec3(1, 0, 0));
    if (ty0 > 0 && ( tmin<0 || ty0 < tmin ) ){ tmin = ty0; }
    tmax = max(tmax, ty0);
  }
  {
    float ty1 = squareRayTest(org, dir, vec3(0, +hlen, 0), hlen, vec3(0, +1, 0), vec3(0, 0, 1), vec3(1, 0, 0));
    if (ty1 > 0 && ( tmin<0 || ty1 < tmin ) ){ tmin = ty1; }
    tmax = max(tmax, ty1);
  }
}

void rayFromCamera(
  out vec3 org,
  out vec3 dir )
{
  vec2 pos = vec2(
    (gl_FragCoord.x*2.0 - resolution.x) / resolution.x,
    (gl_FragCoord.y*2.0 - resolution.y) / resolution.y );

  vec4 p0 = mMVPinv * vec4(pos.x,pos.y,-1,1);
  vec4 p1 = mMVPinv * vec4(pos.x,pos.y,+1,1);
  vec3 camPos0 = p0.xyz/p0.w;
  vec3 camPos1 = p1.xyz/p1.w;
  org = camPos0;
  dir = normalize(camPos1-camPos0);
}

void main() {
  vec3 org,dir;
  rayFromCamera(org,dir);

  float tmin, tmax;
  cubeRayTest(tmin,tmax,org,dir,0.5);

  float len = 0;
  if( tmax > 0 ){
    if( tmin > 0 ){
      len = tmax-tmin;
    }
    else {
      len = tmax;
    }
  }
  float val = 1-exp(-len*5);
  FragColor = vec4(val,val,val,0);
}
