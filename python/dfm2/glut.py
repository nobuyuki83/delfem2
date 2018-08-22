from OpenGL.GL import *
from OpenGL.GLUT import *
import dfm2

class WindowManager:
  def __init__(self,view_height):
    self.camera = dfm2.gl.Camera(view_height)
    self.modifier = 0
    self.mouse_x = 0.0
    self.mouse_y = 0.0

  def special(self,key,x,y):
    if key == int(GLUT_KEY_PAGE_UP):
      self.camera.scale *= 1.03
    elif key == int(GLUT_KEY_PAGE_DOWN):
      self.camera.scale *= (1.0 / 1.03)
    glutPostRedisplay()

  def mouse(self,button, state, x, y):
    self.modifier = glutGetModifiers()
    self.mouse_x, self.mouse_y = dfm2.gl.mouse_screen_pos(x, y)
    glutPostRedisplay()

  def motion(self,x, y):  
    if self.modifier == 1: # shift
      self.mouse_x, self.mouse_y, self.camera.quat, self.camera.trans = dfm2.gl.motion_rot(
          x, y, self.mouse_x, self.mouse_y, self.camera.quat,
          self.camera.trans, self.camera.view_height)
    if self.modifier == 2 or self.modifier == 4: # ctrl or cmnd+alt
      self.mouse_x, self.mouse_y, self.camera.quat, self.camera.trans = dfm2.gl.motion_trans(
          x, y, self.mouse_x, self.mouse_y, self.camera.quat,
          self.camera.trans, self.camera.view_height)
    glutPostRedisplay()


