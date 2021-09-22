#ifndef GLUT_FUNCS_H
#define GLUT_FUNCS_H

#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>  // for the format

#include "delfem2/cam_modelview.h"
#include "delfem2/cam_projection.h"

#if defined(__APPLE__) && defined(__MACH__)
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

// end of header
// --------------------------------------------

static void RenderBitmapString(float x, float y, void *font,const char *string)
{
  ::glRasterPos2f(x, y);
  for(unsigned int i=0;string[i]!='\0';i++){
    ::glutBitmapCharacter(font, string[i]);
  }
}


static void DrawMessage_LeftTop(const char* msg, int iline=0)
{
  int* font=(int*)GLUT_BITMAP_8_BY_13;
  ////
  glPushAttrib(GL_TRANSFORM_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT | GL_POINT_BIT);
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_TEXTURE_2D);
  ::glDisable(GL_DEPTH_TEST);
  GLint viewport[4]; ::glGetIntegerv(GL_VIEWPORT,viewport);
  const int win_w = viewport[2];
  const int win_h = viewport[3];
  ////
  ::glMatrixMode(GL_PROJECTION);
  ::glPushMatrix();
  ::glLoadIdentity();
  ::gluOrtho2D(0, win_w, 0, win_h);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glLoadIdentity();
  ////
  ::glScalef(1, -1, 1);
  ::glTranslatef(0, -(float)win_h, 0);
//  ::glColor3d(0.0, 0.0, 0.0);
  RenderBitmapString(10.0f,18.0f+iline*12.0f, (void*)font, msg);
  ////
  ::glMatrixMode(GL_PROJECTION);
  ::glPopMatrix();
  ::glMatrixMode(GL_MODELVIEW);
  ::glPopMatrix();
  ::glPopAttrib();
}



static void ShowFPS(){
  ::glDisable(GL_LIGHTING);
//  static char s_fps[32];
  static std::string s_fps;
  {
    static int frame, timebase;
    int time;
    frame++;
    time=glutGet(GLUT_ELAPSED_TIME);
    if (time - timebase > 500) {
      std::stringstream ss;
      ss << "FPS:";
      ss << std::fixed << std::setprecision(2);
      ss << frame*1000.0/(time-timebase);
      s_fps = ss.str();
      timebase = time;
      frame = 0;
    }
  }
//  DrawMessage_LeftTop(s_fps);
  DrawMessage_LeftTop(s_fps.c_str());
}



/**
 * @brief class for each GLUT window
 */
class CNav3D_GLUT
{
public:  
  void glutMouse(int button, int state, int x, int y){
    GLint viewport[4];
    ::glGetIntegerv(GL_VIEWPORT,viewport);
    imodifier = glutGetModifiers();
    if( state == GLUT_UP ){ this->ibutton = -1; }
    else{                   this->ibutton = button; }
    const int win_w = viewport[2];
    const int win_h = viewport[3];
    mouse_x = (2.0*x-win_w)/win_w;
    mouse_y = (win_h-2.0*y)/win_h;
    if( state == GLUT_DOWN ){
      mouse_x_down = mouse_x;
      mouse_y_down = mouse_y;
    }
    ::glutPostRedisplay();
  }
  void glutMotion( int x, int y ){
    GLint viewport[4];
    ::glGetIntegerv(GL_VIEWPORT,viewport);
    const int win_w = viewport[2];
    const int win_h = viewport[3];
    const double mov_end_x = (2.0*x-win_w)/win_w;
    const double mov_end_y = (win_h-2.0*y)/win_h;
    dx = mov_end_x - mouse_x;
    dy = mov_end_y - mouse_y;
    {
      if(      imodifier == GLUT_ACTIVE_ALT   ){
        modelview.Rot_Camera(dx, dy);
      }
      else if( imodifier == GLUT_ACTIVE_SHIFT ){
        double s0 = projection.view_height / scale;
        trans[0] += s0 * dx;
        trans[1] += s0 * dy;
      }
    }
    mouse_x = mov_end_x;
    mouse_y = mov_end_y;
    ::glutPostRedisplay();
  }
  void glutSpecial(int Key, [[maybe_unused]] int x, [[maybe_unused]] int y)
  {
    switch(Key)
    {
      case GLUT_KEY_PAGE_UP:
        scale *= 1.03;
        break;
      case GLUT_KEY_PAGE_DOWN:
        scale *= (1.0/1.03);
        break;
      case GLUT_KEY_F1:
        projection.is_pars = !projection.is_pars;
        break;
      case GLUT_KEY_F2:
        projection.fovy *= 1.05;
        break;
      case GLUT_KEY_F3:
        projection.fovy /= 1.05;
        break;
      case GLUT_KEY_LEFT:
        break;
      case GLUT_KEY_RIGHT:
        break;
      case GLUT_KEY_UP:
        break;
      case GLUT_KEY_DOWN:
        break;
      case GLUT_KEY_HOME :
        break;
      case GLUT_KEY_END :
        break;
      default:
        break;
    }
    ::glutPostRedisplay();
  }
  void SetGL_Camera()
  {
    GLint viewport[4]; ::glGetIntegerv(GL_VIEWPORT,viewport);
    const int win_w = viewport[2];
    const int win_h = viewport[3];
    {
      ::glMatrixMode(GL_PROJECTION);
      ::glLoadIdentity();
      double asp = (double)win_w/win_h;
      const delfem2::CMat4f mP = projection.GetMatrix(asp);
      const delfem2::CMat4f mS = delfem2::CMat4f::Scale(scale);
      // const delfem2::CMat4f mZ = delfem2::CMat4f::ScaleXYZ(1,1,-1);
      ::glMultMatrixf((mS * mP.transpose()).data());
    }
    {
      ::glMatrixMode(GL_MODELVIEW);
      ::glLoadIdentity();
      const delfem2::CMat4f mMV = modelview.Mat4ColumnMajor();
      const delfem2::CMat4f mT = delfem2::CMat4f::Translate(trans);
      ::glMultMatrixf((mMV*mT.transpose()).data());
    }
  }
public:
  int iwin;
  int imodifier;
  int ibutton;
  delfem2::Projection_LookOriginFromZplus<double> projection;
  delfem2::ModelView_Trackball<double> modelview;
  double scale = 1.0;
  double trans[3] = {0,0,0};
  double mouse_x, mouse_y;
  double dx;
  double dy;
  double mouse_x_down, mouse_y_down;
};



#endif
