varying vec2 vST;
varying vec3 lightVec; 
varying vec3 eyeVec;

void main()
{
	vec3 vTangent = vec3(0.0);
  vec3 c1 = cross(gl_Normal, vec3(0.0, 0.0, 1.0)); 
  vec3 c2 = cross(gl_Normal, vec3(0.0, 1.0, 0.0)); 
  if(length(c1)>length(c2)){ vTangent = c1;	}
  else{                      vTangent = c2;	}
  vTangent = normalize(vTangent);

  vec3 n = normalize(gl_NormalMatrix * gl_Normal);
  vec3 t = normalize(gl_NormalMatrix * vTangent);
  vec3 b = cross(n, t);

  vec3 vVertex = vec3(gl_ModelViewMatrix * gl_Vertex);
  vec3 tmpVec = gl_LightSource[0].position.xyz - vVertex;

//  tmpVec = gl_NormalMatrix * tmpVec;

  lightVec.x = dot(tmpVec, t);
  lightVec.y = dot(tmpVec, b);
  lightVec.z = dot(tmpVec, n);

  tmpVec = -vVertex;
  eyeVec.x = dot(tmpVec, t);
  eyeVec.y = dot(tmpVec, b);
  eyeVec.z = dot(tmpVec, n);

	vST = gl_MultiTexCoord0.st;
  gl_Position = ftransform();
}