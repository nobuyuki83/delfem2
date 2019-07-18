#ifndef FUNCS_GLFW_H
#define FUNCS_GLFW_H

#include "gl_camera.h"

// class for each GLUT window
class CWindowManager_GLFW
{
public:
  void Mouse(GLFWwindow *window, int button, int action, int mods){
    GLint viewport[4];
    ::glGetIntegerv(GL_VIEWPORT,viewport);
//    int win_w, win_h;  glfwGetFramebufferSize(window, &win_w, &win_h);
    int win_w = 800;//viewport[2];
    int win_h = 500;//viewport[3];
    imodifier = mods;
    double x, y;  glfwGetCursorPos (window, &x,&y);
//    std::cout << " pos: " << x << " " << y << " " << win_w << " " << win_h << std::endl;
    mouse_x = (2.0*x-win_w)/win_w;
    mouse_y = (win_h-2.0*y)/win_h;
    if( action == 0 ){
      ibutton = -1;
    }
    else if( action == 1 ){ // mouse down
      ibutton = button;
      mouse_x_down = mouse_x;
      mouse_y_down = mouse_y;
    }
  }
  void Motion(GLFWwindow *window, double x, double y){

    int win_w, win_h;  glfwGetFramebufferSize(window, &win_w, &win_h);
    const double mov_end_x = (2.0*x-win_w)/win_w;
    const double mov_end_y = (win_h-2.0*y)/win_h;
    dx = mov_end_x - mouse_x;
    dy = mov_end_y - mouse_y;
    if( ibutton == -1 ){
    }
    else{
      if(      imodifier == 4  ){
        camera.Rot_Camera(dx, dy);
      }
      else if( imodifier == 1 ){
        camera.Pan_Camera(dx, dy);
      }
    }
    mouse_x = mov_end_x;
    mouse_y = mov_end_y;
  }
  void Key(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    if( key ==  GLFW_KEY_PAGE_UP && action == GLFW_PRESS ){ camera.Scale(1.03); }
    if( key ==  GLFW_KEY_PAGE_DOWN && action == GLFW_PRESS ){ camera.Scale(1.0/1.03); }
  }
  void SetGL_Camera(GLFWwindow *window){
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    camera.SetGL_Camera(width,height);
  }
public:
//  int iwin;
  int imodifier;
  int ibutton;
  CCamera camera;
  double mouse_x, mouse_y;
  double dx, dy;
  double mouse_x_down, mouse_y_down;
};

#endif /* utility_glfw_h */
