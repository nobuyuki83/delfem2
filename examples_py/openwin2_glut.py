from OpenGL.GL import *
from OpenGL.GLUT import *

import sys
sys.path.append("../module_py")
import dfm2
import dfm2.glut

def draw_func():
  glEnable(GL_LIGHTING)
  glutSolidTeapot(0.5)

win = dfm2.glut.WindowGLUT(1.0,winsize=(400,300))
dfm2.setSomeLighting()
win.draw_loop(draw_func)

