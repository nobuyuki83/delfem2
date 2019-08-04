####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import OpenGL.GL as gl
import glfw

import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

msh = None
nav = None

def mouseButtonCB(win_glfw,btn,action,mods):
  nav.mouse(win_glfw, btn, action, mods)

def mouseMoveCB(win_glfw,x,y):
  nav.motion(win_glfw, x, y)

def keyFunCB(win_glfw,key,scancode,action,mods):
  nav.keyinput(win_glfw, key, scancode, action, mods)

def render():
  gl.glColor3d(1, 0, 0)
  gl.glEnable(gl.GL_LIGHTING)
  msh.draw()

def main():
  global msh, nav
  msh = dfm2.Mesh()
  msh.read("../test_inputs/bunny_2k.ply")
  msh.scale_xyz(0.02)

  nav = dfm2.gl.glfw.NavigationGLFW(1.0)

  glfw.init()
  win_glfw = glfw.create_window(640, 480, 'Hello World', None, None)
  glfw.make_context_current(win_glfw)

  dfm2.gl.setSomeLighting()
  gl.glEnable(gl.GL_DEPTH_TEST)

  glfw.set_mouse_button_callback(win_glfw, mouseButtonCB)
  glfw.set_cursor_pos_callback(win_glfw,  mouseMoveCB)
  glfw.set_key_callback(win_glfw, keyFunCB)

  while not glfw.window_should_close(win_glfw):
    gl.glClearColor(1, 1, 1, 1)
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
    gl.glPolygonOffset(1.1, 4.0)

    nav.camera.set_gl_camera()
    render()
    glfw.swap_buffers(win_glfw)
    glfw.poll_events()
    if nav.isClose:
      break
  glfw.destroy_window(win_glfw)
  glfw.terminate()
  print("closed")

if __name__ == "__main__":
  main()