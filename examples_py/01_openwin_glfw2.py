####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import OpenGL.GL as gl

import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw


def draw_func():
  gl.glEnable(gl.GL_LIGHTING)
  msh.draw()

msh = dfm2.Mesh()
msh.read("../test_inputs/bunny_2k.ply")
msh.scale_xyz(0.03)

win = dfm2.gl.glfw.WindowGLFW(1.0,winsize=(400,300))
win.list_func_draw.append(draw_func)
dfm2.gl.setSomeLighting()
win.draw_loop()
