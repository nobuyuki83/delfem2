#include <iostream>
#include <vector>
#include <set>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/glut_funcs.h"
#include "delfem2/gl_funcs.h"

#include "delfem2/funcs.h"
#include "delfem2/mshio.h"
#include "delfem2/../../external/io_fbx.h"

#define STB_IMAGE_IMPLEMENTATION
#include "delfem2/../../external/stb/stb_image.h"

CGlutWindowManager window;
CTexManager tex_manager;
CRigMsh rigmsh;

void myGlutDisplay(void)
{
  //  ::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(0.5f, 0.8f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  window.SetGL_Camera();
  
  
  rigmsh.Draw(tex_manager);
  
  DrawAxis(100);
  
  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}


void myGlutIdle()
{
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  ::glViewport(0,0,w,h);
  ::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  window.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  window.glutMotion(x, y);
  if( window.imodifier != 0 ) return;
  ////
  rigmsh.Drag(window.mouse_x,window.mouse_y,window.dx,window.dy);
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  window.glutMouse(button, state, x, y);
  /////
  window.SetGL_Camera();
  if( state == GLUT_DOWN ){
    rigmsh.Pick(window.mouse_x,window.mouse_y);
  }
  else{
    rigmsh.ielem_selected = -1;
  }
  /////
  ::glutPostRedisplay();
}

void Read(const std::string& path_fbx){
  Read_FBX(path_fbx,rigmsh);
  for(int ib=0;ib<rigmsh.aBone.size();++ib){
    std::cout << ib << " " << rigmsh.aBone[ib].name << std::endl;
  }
  rigmsh.FixTexPath(path_fbx);
  rigmsh.PrintInfo();
  { // go all the textures
    std::vector<std::string> aPath = rigmsh.GetArrayTexPath();
    tex_manager.Clear();
    for(int ipath=0;ipath<aPath.size();++ipath){
      const std::string path = aPath[ipath];
      int width, height, bpp;
      unsigned char* pixels; pixels = stbi_load(path.c_str(), &width, &height, &bpp, 0);
      if( pixels == 0){ continue; }
      std::cout << "Texture Loaded: " << path << std::endl;
      stbi__vertical_flip(pixels, width, height, bpp);
      tex_manager.AddTexture(pixels, path, width, height, bpp);
      stbi_image_free(pixels);
    }
  }
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'm':
    {
      break;
    }
    case 'b':
    {
      rigmsh.is_draw_bone = !rigmsh.is_draw_bone;
      break;
    }
    case 'W':
    {
      int nbone = rigmsh.aBone.size();
      rigmsh.ibone_selected = (rigmsh.ibone_selected+1)%(nbone+1);
      break;
    }
    case '1':
      break;
    case '2':
      break;
    case '3':
      break;
    case '4':
      break;
    case 'i': // one iteration
      break;
    case 'd': // change draw mode
      break;
    case 'f': //
      break;
    case 's': //
      break;
  }
  ::glutPostRedisplay();
}


int main(int argc,char* argv[])
{
  glutInit(&argc, argv);
  
  // Initialize GLUT window 3D
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_STENCIL);
  glutCreateWindow("3D View");
  glutDisplayFunc(myGlutDisplay);
  glutIdleFunc(myGlutIdle);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  
  ////////////////////////
  
  window.camera.view_height = 100.0;
  window.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  
  std::cout << argc << std::endl;
  if( argc < 2 ){
    std::cout << "require path to the FBX in the argument" << std::endl;
    return 0;
  }
  
  std::cout << "open fbx. path: " << argv[1] << std::endl;
  if( !isFileExists(argv[1]) ){
    std::cout << "file doesn't exist" << std::endl;
    return 0;
  }
  
  Read(argv[1]);

  
  glutMainLoop();
  return 0;
}
