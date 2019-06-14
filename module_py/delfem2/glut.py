####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


from OpenGL.GL import *
from OpenGL.GLUT import *
from .gl import *

def draw_text(x, y, font, text, color):
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glMatrixMode(GL_MODELVIEW)
  glLoadIdentity()
  glColor3f(color[0], color[1], color[2])
  glRasterPos2f(x, y)
  for ch in text:
    glutBitmapCharacter(font, ctypes.c_int(ord(ch)))

def draw_sphere(pos, rad, color):
  if pos is None:
    return
  glColor3d(color[0], color[1], color[2])
  glTranslatef(+pos[0], +pos[1], +pos[2])
  glutSolidSphere(rad, 32, 32)
  glTranslatef(-pos[0], -pos[1], -pos[2])


class WindowManagerGLUT:
  def __init__(self,view_height):
    self.camera = Camera(view_height)
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
    viewport = glGetIntegerv(GL_VIEWPORT)
    (win_w,win_h) = viewport[2:4]
    self.mouse_x = (2.0 * x - win_w) / win_w
    self.mouse_y = (win_h - 2.0 * y) / win_h
    self.modifier = glutGetModifiers()
    glutPostRedisplay()

  def motion(self,x, y):  
    viewport = glGetIntegerv(GL_VIEWPORT)
    if self.modifier == 2 or self.modifier == 4:  # ctrl or cmnd+alt
      self.mouse_x, self.mouse_y = self.camera.rotation(
          x, y, self.mouse_x, self.mouse_y,
          viewport[2],viewport[3])
    ####
    if self.modifier == 1:  # shift
      self.mouse_x, self.mouse_y = self.camera.translation(
          x, y, self.mouse_x, self.mouse_y,
          viewport[2],viewport[3])
    glutPostRedisplay()

class WindowGLUT:
  def __init__(self,view_height,winsize=(400,300)):
    self.wm = WindowManagerGLUT(view_height)
    self.draw_func = None
    glutInit()
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)  # zBuffer
    glutInitWindowSize(winsize[0], winsize[1])
    glutInitWindowPosition(100, 100)
    glutCreateWindow("Visualizzatore_2.0")

  def keyboard(self, bkey, x, y):
    key = bkey.decode("utf-8")
    if key == 'q':
      exit()
    glutPostRedisplay()      

  def display(self):
    glClearColor(0.3, 0.5, 0.8, 1.0)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glEnable(GL_DEPTH_TEST)  
    self.wm.camera.set_gl_camera()
    self.draw_func()
    glutSwapBuffers()  

  def reshape(self, width, height):
    glViewport(0, 0, width, height)

  def idle(self):
    glutPostRedisplay()

  def special(self,key, x, y):  
    self.wm.special(key,x,y)
    glutPostRedisplay()    

  def mouse(self,button, state, x, y):
    self.wm.mouse(button,state,x,y)
    glutPostRedisplay()    

  def motion(self,x, y):
    self.wm.motion(x,y)    
    glutPostRedisplay()  

  def draw_loop(self,draw_func0):  

    # Register callbacks
    glutReshapeFunc(self.reshape)
    glutDisplayFunc(self.display)
    glutMouseFunc(self.mouse)
    glutMotionFunc(self.motion)
    glutKeyboardFunc(self.keyboard)
    glutSpecialFunc(self.special)
    glutIdleFunc(self.idle)

    self.draw_func = draw_func0
    glutMainLoop()


