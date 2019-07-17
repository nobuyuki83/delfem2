####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import OpenGL.GL as gl
import OpenGL.GLUT as glut

from ._gl import Camera, screenUnProjection, screenUnProjectionDirection, AxisXYZ

from .c_gl import setSomeLighting



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
    self.ratio_x1 = 0.0
    self.ratio_y1 = 0.0
    self.ratio_x0 = 0.0
    self.ratio_y0 = 0.0

  def special(self,key,x,y):
    if key == int(GLUT_KEY_PAGE_UP):
      self.camera.scale *= 1.03
    elif key == int(GLUT_KEY_PAGE_DOWN):
      self.camera.scale *= (1.0 / 1.03)
    glutPostRedisplay()

  def mouse(self,button, state, x, y):
    viewport = gl.glGetIntegerv(gl.GL_VIEWPORT)
    (win_w,win_h) = viewport[2:4]
    self.ratio_x1 = (2.0 * x - win_w) / win_w
    self.ratio_y1 = (win_h - 2.0 * y) / win_h
    self.modifier = glut.glutGetModifiers()
    glut.glutPostRedisplay()

  def motion(self,x, y):
    self.ratio_x0,self.ratio_y0 = self.ratio_x1,self.ratio_y1
    viewport = gl.glGetIntegerv(gl.GL_VIEWPORT)
    (win_w,win_h) = viewport[2:4]
    self.ratio_x1 = (2.0 * x - win_w) / win_w
    self.ratio_y1 = (win_h - 2.0 * y) / win_h
    if self.modifier == 2 or self.modifier == 4:  # ctrl or cmnd+alt
      self.camera.rotation(
          self.ratio_x1, self.ratio_y1, self.ratio_x0, self.ratio_y0)
    ####
    if self.modifier == 1:  # shift
      self.camera.translation(
        self.ratio_x1, self.ratio_y1, self.ratio_x0, self.ratio_y0)
    glut.glutPostRedisplay()

class WindowGLUT:
  def __init__(self,view_height,winsize=(400,300)):
    self.wm = WindowManagerGLUT(view_height)
    self.draw_func = None
    glut.glutInit()
    glut.glutInitDisplayMode(glut.GLUT_DOUBLE | glut.GLUT_RGB | glut.GLUT_DEPTH)  # zBuffer
    glut.glutInitWindowSize(winsize[0], winsize[1])
    glut.glutInitWindowPosition(100, 100)
    glut.glutCreateWindow("Visualizzatore_2.0")

  def keyboard(self, bkey, x, y):
    key = bkey.decode("utf-8")
    if key == 'q':
      exit()
    glutPostRedisplay()      

  def display(self):
    gl.glClearColor(0.3, 0.5, 0.8, 1.0)
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    gl.glEnable(gl.GL_DEPTH_TEST)
    self.wm.camera.set_gl_camera()
    self.draw_func()
    glut.glutSwapBuffers()

  def reshape(self, width, height):
    gl.glViewport(0, 0, width, height)

  def idle(self):
    glut.glutPostRedisplay()

  def special(self,key, x, y):  
    self.wm.special(key,x,y)
    glut.glutPostRedisplay()

  def mouse(self,button, state, x, y):
    self.wm.mouse(button,state,x,y)
    glut.glutPostRedisplay()

  def motion(self,x, y):
    self.wm.motion(x,y)    
    glut.glutPostRedisplay()

  def draw_loop(self,draw_func0):  

    # Register callbacks
    glut.glutReshapeFunc(self.reshape)
    glut.glutDisplayFunc(self.display)
    glut.glutMouseFunc(self.mouse)
    glut.glutMotionFunc(self.motion)
    glut.glutKeyboardFunc(self.keyboard)
    glut.glutSpecialFunc(self.special)
    glut.glutIdleFunc(self.idle)

    self.draw_func = draw_func0
    glut.glutMainLoop()


