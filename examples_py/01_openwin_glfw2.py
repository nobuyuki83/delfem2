####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw


def draw_func():
  glEnable(GL_LIGHTING)
  msh.draw()

msh = dfm2.Mesh()
msh.read("../test_inputs/bunny_2k.ply")
msh.scale_xyz(0.03)

win = dfm2.glfw.WindowGLFW(1.0,winsize=(400,300))
win.list_func_draw.append(draw_func)
dfm2.setSomeLighting()
win.draw_loop()
