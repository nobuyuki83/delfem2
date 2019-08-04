####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


from OpenGL.GL import *
from OpenGL.GLUT import *

import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glut

def draw_func():
  glEnable(GL_LIGHTING)
  glutSolidTeapot(0.5)

win = dfm2.gl.glut.WindowGLUT(1.0,winsize=(400,300))
dfm2.gl.setSomeLighting()
win.draw_loop(draw_func)

