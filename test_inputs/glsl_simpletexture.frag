#version 120

vec3 ADSLightModel(
  vec3 myNormal, vec3 myPosition, 
  gl_MaterialParameters myMaterial,
  int indexLight)
{
  vec3 norm = normalize(myNormal);
  vec3 light_pos0 = gl_LightSource[indexLight].position.xyz;
  vec3 lightv = normalize( light_pos0 - myPosition ); 
  vec3 ambient = myMaterial.ambient.xyz * gl_LightSource[indexLight].ambient.xyz;
  vec3 specular = vec3(0.0, 0.0, 0.0);
  {
    vec3 viewv  = normalize( -myPosition ); 
    vec3 refl   = reflect(-lightv,norm); 
    if( dot(lightv,viewv) > 0.0 ){
      float sh0 = myMaterial.shininess;
      specular = pow( max(dot(refl,viewv), 0.0), sh0) * myMaterial.specular.xyz;
    }   
  }
  vec3 diffuse = max(dot(lightv,norm),0.0) * myMaterial.diffuse.xyz;
  return clamp( specular+diffuse+ambient, 0.0, 1.0);
}

varying vec2 vST;
varying vec3 lightVec; 
varying vec3 eyeVec;
uniform sampler2D myTexColor;
uniform sampler2D myTexNormal;
void main()
{
//	sampler2D sl0 = 0;
  vec3 c0 = texture2D( myTexColor,  vST ).rgb;
  vec3 n0 = texture2D( myTexNormal, vST ).rgb;
  n0 = normalize(n0 * 2.0 - 1.0);

//  vec3 light_pos0 = gl_LightSource[indexLight].position.xyz;
  vec3 lightv = vec3(0., 0., 1.0);
//  vec3 lightv = normalize( light_pos0 ); 
//  gl_FragColor = vec4(c0,1);
//	vec3 c0 = vec3(1.0, 1.0, 0.0);
  vec3 c1 = max(dot(n0,normalize(lightVec)),0.0) * c0;
//  vec3 c1 = mix(c0, n0, 0.5);
  gl_FragColor = vec4(c1, 1.0);
}