varying vec3 LightIntensity; 

vec3 ADSLightModel(
  vec3 myNormal, vec3 myPosition, 
  gl_MaterialParameters myMaterial,
  int indexLight)
{    
  vec3 norm = normalize(myNormal);
  vec3 light_pos0 = gl_LightSource[indexLight].position.xyz;
  vec3 lightv = normalize( light_pos0 - myPosition ); 
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
  return clamp( specular+diffuse, 0.0, 1.0);
}

void main(){ 
  gl_Position = ftransform();
  vec3 myNormal   = vec3(gl_NormalMatrix    * gl_Normal);
  vec3 myPosition = vec3(gl_ModelViewMatrix * gl_Vertex);
  LightIntensity  = ADSLightModel(myNormal,myPosition,gl_FrontMaterial,0);
}
