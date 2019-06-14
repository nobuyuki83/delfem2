####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


from OpenGL.GL import *
from OpenGL.GLUT import *

import sys
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glut

def draw_func():
  glEnable(GL_LIGHTING)
  glutSolidTeapot(0.5)

win = dfm2.glut.WindowGLUT(1.0,winsize=(400,300))
dfm2.setSomeLighting()
win.draw_loop(draw_func)

