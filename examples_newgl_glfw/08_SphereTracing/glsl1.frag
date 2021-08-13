#version 330 core

uniform vec2 resolution;
uniform mat4 mMVPinv;
uniform mat4 mMV;

out vec4 FragColor;

vec3 trans(vec3 p){
  return mod(p+1.5, 3.) - 1.5;
}

float distanceFunction(vec3 pos0){
//  return length(trans(pos)) - 1.0;
  vec3 pos1 = trans(pos0);
  float x0 = max(abs(pos1.x),abs(pos1.y));
  return max(x0,abs(pos1.z))-0.5;
}

vec3 getNormal(vec3 p)
{
  const float d = 0.0001;
  float nx0 = distanceFunction(p+vec3(d, 0.0, 0.0))-distanceFunction(p+vec3(-d, 0.0, 0.0));
  float ny0 = distanceFunction(p+vec3(0.0, d, 0.0))-distanceFunction(p+vec3(0.0, -d, 0.0));
  float nz0 = distanceFunction(p+vec3(0.0, 0.0, d))-distanceFunction(p+vec3(0.0, 0.0, -d));
  return normalize(vec3(nx0, ny0, nz0));
}

bool sphere_tracing(
  inout vec3 normal,
  vec3 camPos1,
  vec3 camDir1)
{ // sphere tracing
  vec3 posOnRay = camPos1;
  float distance;
  for (int i=0; i<64; ++i){
    distance = distanceFunction(posOnRay);
    posOnRay += distance*camDir1;
  }
  normal = getNormal(posOnRay);
  return abs(distance) < 0.001;
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

  vec3 normal;
  bool is_hit = sphere_tracing(normal, org, dir);
  normal = (mMV * vec4(normal,0)).xyz;

  if( is_hit ){
    FragColor = vec4(normal*0.5+0.5, 1.0);
  }
  else{
    FragColor = vec4(1.0);
  }
}
