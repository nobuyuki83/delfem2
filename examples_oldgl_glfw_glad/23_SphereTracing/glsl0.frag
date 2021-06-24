uniform vec2 resolution;

struct Ray{
  vec3 origin;
  vec3 direction;
};

struct Sphere{
  float radius;
  vec3  position;
  vec3  color;
};

struct Intersection{
  bool hit;
  vec3 hitPoint;
  vec3 normal;
  vec3 color;
};

Intersection intersectSphere(Ray R, Sphere S){
  Intersection i;
  vec3  a = R.origin - S.position;
  float b = dot(a, R.direction);
  float c = dot(a, a) - (S.radius * S.radius);
  float d = b * b - c;
  if(d > 0.0){
    float t = -b - sqrt(d);
    if(t > 0.0){
      i.hit = true;
      i.hitPoint = R.origin + R.direction * t;
      i.normal = normalize(i.hitPoint - S.position);
      float d = clamp(dot(normalize(vec3(1.0)), i.normal), 0.1, 1.0);
      i.color = S.color * d;
      return i;
    }
  }
  i.hit = false;
  i.hitPoint = vec3(0.0);
  i.normal = vec3(0.0);
  i.color = vec3(0.0);
  return i;
}

void main()
{
  vec2 p = (gl_FragCoord.xy * 2.0 - resolution) / min(resolution.x, resolution.y);
  
  Ray ray;
  ray.origin = vec3(p.x, p.y, 5.0);
  ray.direction = vec3(0, 0, -1.0);
  
  Sphere sphere;
  sphere.radius = 1.0;
  sphere.position = vec3(0.0);
  sphere.color = vec3(1.0);
  
  // hit check
  Intersection i = intersectSphere(ray, sphere);
  
  if(i.hit){
    gl_FragColor = vec4( i.color, 1.0 );
  }
  else{
    gl_FragColor = vec4(vec3(0,0,0), 1.0);
  }
  
}
