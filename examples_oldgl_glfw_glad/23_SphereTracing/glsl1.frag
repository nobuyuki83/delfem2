uniform vec2 resolution;
uniform float time;

vec3 trans(vec3 p){
  return mod(p, 4.0) - 2.0;
}

float distanceFunction(vec3 pos){
  return length(trans(pos)) - 1.0;
}
/*
float distanceFunction(vec3 pos)
{
  return length(pos) - 1.0;
}
 */

vec3 getNormal(vec3 p)
{
  const float d = 0.0001;
  return
  normalize
  (
  vec3
  (
  distanceFunction(p+vec3(d,0.0,0.0))-distanceFunction(p+vec3(-d,0.0,0.0)),
  distanceFunction(p+vec3(0.0,d,0.0))-distanceFunction(p+vec3(0.0,-d,0.0)),
  distanceFunction(p+vec3(0.0,0.0,d))-distanceFunction(p+vec3(0.0,0.0,-d))
  )
  );
}

void main() {
  vec2 pos = (gl_FragCoord.xy*2.0 -resolution) / resolution.y;
  
  vec3 camPos = vec3(0.0, 0.0, 3.0-time);
  vec3 camDir = vec3(0.1*sin(time), 0.2*cos(time), -1.0);
  vec3 camUp = vec3(0.0, 1.0, 0.0);
  vec3 camSide = cross(camDir, camUp);
  float focus = 1.8;
  
  vec3 rayDir = normalize(camSide*pos.x + camUp*pos.y + camDir*focus);
  
  float t = 0.0, d;
  vec3 posOnRay = camPos;
  
  for(int i=0; i<64; ++i)
  {
    d = distanceFunction(posOnRay);
    t += d;
    posOnRay = camPos + t*rayDir;
  }
  
  vec3 normal = getNormal(posOnRay);
  if(abs(d) < 0.001)
  {
    gl_FragColor = vec4(normal, 1.0);
  }else
  {
    gl_FragColor = vec4(0.0);
  }
}
