####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import OpenGL.GL as gl
import numpy

import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def main():
  buf_face = None
  buf_edge = None

  def draw_func():
    gl.glEnable(gl.GL_LIGHTING)
    buf_face.draw()
    gl.glDisable(gl.GL_LIGHTING)
    gl.glColor3d(0,0,0)
    buf_edge.draw()

  msh = dfm2.Mesh()
  msh.read("../test_inputs/bunny_2k.ply")
  msh.scale_xyz(0.03)
  msh_edge = msh.mesh_edge()
  ####
  win = dfm2.gl.glfw.WindowGLFW(1.0,winsize=(400,300))
  buf_face = dfm2.gl.GLBufferMesh(msh,is_normal=True)
  buf_edge = dfm2.gl.GLBufferMesh(msh_edge)
  win.list_func_draw.append(draw_func)
  dfm2.gl.setSomeLighting()
  win.draw_loop()

if __name__ == "__main__":
  main()