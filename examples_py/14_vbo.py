####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import OpenGL.GL as gl
import numpy

import sys
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw


msh = None
vbo_xy = None
ebo_tri = None

def vbo_array(aXY:numpy.ndarray) -> int:
  assert aXY.dtype == numpy.float64
  vbo = gl.glGenBuffers(1)
  gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)
  gl.glBufferData(gl.GL_ARRAY_BUFFER,
                  aXY.size * 8,
                  (gl.ctypes.c_double * aXY.size)(*aXY.flatten()),
                  gl.GL_STATIC_DRAW)
  gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
  return vbo

def ebo_array(aElm:numpy.ndarray) -> int:
  ebo = gl.glGenBuffers(1)
  gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, ebo)
  gl.glBufferData(gl.GL_ELEMENT_ARRAY_BUFFER,
                  aElm.size*4,
                  (gl.ctypes.c_uint32 * aElm.size)(*aElm.flatten()),
                  gl.GL_STATIC_DRAW)
  gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, 0)
  print(ebo,aElm,type(aElm),aElm.dtype)
  return ebo


def draw(msh,vbo_xy,ebo_tri):
  gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
  gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo_xy)
  gl.glVertexPointer(msh.np_pos.shape[1], gl.GL_DOUBLE, 0, None)
  gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, ebo_tri)
  gl.glDrawElements(gl.GL_TRIANGLES, msh.np_elm.size, gl.GL_UNSIGNED_INT, None)
  gl.glDisableClientState(gl.GL_VERTEX_ARRAY)


def draw_func():
  gl.glEnable(gl.GL_LIGHTING)
  draw(msh,vbo_xy,ebo_tri)

def main():
  global msh, vbo_xy, ebo_tri
  msh = dfm2.Mesh()
  msh.read("../test_inputs/bunny_2k.ply")
  msh.scale_xyz(0.03)

  win = dfm2.glfw.WindowGLFW(1.0,winsize=(400,300))
  vbo_xy = vbo_array(msh.np_pos)
  ebo_tri = ebo_array(msh.np_elm)

  win.list_func_draw.append(draw_func)
  dfm2.setSomeLighting()
  win.draw_loop()

if __name__ == "__main__":
  main()