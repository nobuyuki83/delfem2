####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import OpenGL.GL as gl

import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

from PyDelFEM2.gl.c_gl import cppDrawSphere_Edge, cppDrawTorus_Edge

def main():
  with dfm2.gl.glfw.WindowGLFW(1.0,winsize=(400,300)) as win0:
    gl.glDisable(gl.GL_LIGHTING)
    gl.glColor3d(0,0,0)
    win0.list_func_draw.append(lambda: cppDrawSphere_Edge(0.5))
    win0.draw_loop()

  with dfm2.gl.glfw.WindowGLFW(1.0,winsize=(400,300)) as win0:
    gl.glDisable(gl.GL_LIGHTING)
    gl.glColor3d(0,0,0)
    win0.list_func_draw.append(lambda: cppDrawTorus_Edge(0.5, 0.3))
    win0.draw_loop()


if __name__ == "__main__":
  main()
