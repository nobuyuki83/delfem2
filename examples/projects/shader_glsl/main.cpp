

#include <iostream>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include "delfem2/funcs_glut.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/funcs_glew.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
CGlutWindowManager win;
int id_shader_program = 0;

/////////

std::string LoadFile
(const std::string& fname)
{
  std::ifstream inputFile1(fname.c_str());
  std::istreambuf_iterator<char> vdataBegin(inputFile1);
  std::istreambuf_iterator<char> vdataEnd;
  return std::string(vdataBegin,vdataEnd);
}

void setShaderProgram(int isp){
  std::string glslVert, glslFrag;
  glUseProgram(0);
  glDeleteProgram(id_shader_program);
  if( isp == 0 ){
    glslVert = LoadFile("../test_inputs/glsl_adsphong.vert");
    glslFrag = LoadFile("../test_inputs/glsl_adsphong.frag");
  }
  else if( isp == 1 ){
    glslVert = LoadFile("../test_inputs/glsl_adsgouraud.vert");
    glslFrag = LoadFile("../test_inputs/glsl_adsgouraud.frag");
  }
  else if( isp == 2 ){
    glslVert = LoadFile("../test_inputs/glsl_toon.vert");
    glslFrag = LoadFile("../test_inputs/glsl_toon.frag");
  }
  else if( isp == 3 ){
    glslVert = LoadFile("../test_inputs/glsl_simpletexture.vert");
    glslFrag = LoadFile("../test_inputs/glsl_simpletexture.frag");
  }
  id_shader_program = setUpGLSL(glslVert, glslFrag);
  ////
  
  glUseProgram(id_shader_program);
  {
    GLint texLoc = glGetUniformLocation(id_shader_program, "myTexColor");
    glUniform1i(texLoc, 0); // GL_TEXTURE0
  }
  
  {
    GLint texLoc = glGetUniformLocation(id_shader_program, "myTexNormal");
    glUniform1i(texLoc, 1); // GL_TEXTURE1
  }
}

/////////////////////////////////////////////////////////

void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClearStencil(0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);
  
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
  ::glEnable(GL_LIGHTING);

  win.SetGL_Camera();
  
  {
    float lightPosition[4] = { 0.0, 0.0, 5.0, 1.0 };
    float light_ambient[4] = { 0.3, 0.3, 0.3, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  }
  {
    float kd[4] = {1.0, 0.0, 0.0, 1.0};
    float shininess = 100.0;
    float ks[4] = {1.0, 1.0, 1.0, 1.0};
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,kd);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,kd);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS,shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, ks);
  }
  
  glUseProgram(id_shader_program);
  ::glutSolidTeapot(1.0);
  
  glUseProgram(0);
  ::glDisable(GL_LIGHTING);
  glColor3d(0,0,0);
  ShowFPS();
  
  ::glutSwapBuffers();
}

void myGlutIdle(){
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  ::glViewport(0,0,w,h);
	::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  win.glutSpecial(Key, x, y);
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  win.glutMotion(x, y);
	::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  win.glutMouse(button, state, x, y);
	::glutPostRedisplay();
}


void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case ' ':
    {
      static int ishader = 0;
      ishader = (ishader+1)%4;
      setShaderProgram(ishader);
      break;
    }
      
  }
  
  ::glutPostRedisplay();
}



int main(int argc,char* argv[])
{
	// Initialize GLUT window 3D
  glutInit(&argc, argv);
  glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
 	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_STENCIL);
  glutCreateWindow("3D View");
	glutDisplayFunc(myGlutDisplay);
	glutIdleFunc(myGlutIdle);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutSpecialFunc(myGlutSpecial);
  glutKeyboardFunc(myGlutKeyboard);

  ////////////////////////
#ifndef __MACH__
  glewInit();
#endif
  
  std::cout<<"Vendor:"<<glGetString(GL_VENDOR)<<'\n';
  std::cout<<"GPU: "<<glGetString(GL_RENDERER)<<'\n';
  std::cout<<"OpenGL ver.: "<<glGetString(GL_VERSION)<<'\n';
  std::cout<<"OpenGL shading ver.: " <<glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;
  
  setShaderProgram(3);
  
  /////////////////////////  
  win.camera.view_height = 2.0;
  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  SFile_TGA tga_color;  LoadTGAFile("../test_inputs/rock_color.tga",  &tga_color);
  SFile_TGA tga_normal; LoadTGAFile("../test_inputs/rock_normal.tga", &tga_normal);
  
  GLuint aIndTex[2];
  ::glGenTextures(2, aIndTex);
  
  // color
  ::glActiveTexture(GL_TEXTURE0);
  ::glBindTexture(GL_TEXTURE_2D, aIndTex[0]);
  ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  ::glTexImage2D(GL_TEXTURE_2D,
                 0, GL_RGBA,
                 tga_color.imageWidth, tga_color.imageHeight,
                 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, tga_color.imageData);
  
  // normal
  ::glActiveTexture(GL_TEXTURE1);
  ::glBindTexture(GL_TEXTURE_2D, aIndTex[1]);
  ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  ::glTexImage2D(GL_TEXTURE_2D,
                 0, GL_RGBA,
                 tga_normal.imageWidth, tga_normal.imageHeight,
                 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, tga_normal.imageData);

  ::glEnable(GL_TEXTURE_2D);
   
  glutMainLoop();
	return 0;
}


