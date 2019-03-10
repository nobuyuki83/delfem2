from OpenGL.GL import *
import glfw

import sys
sys.path.append("../module_py")
import dfm2

msh = None
wmngr_glfw = None

def mouseButtonCB(win_glfw,btn,action,mods):
  wmngr_glfw.mouse(win_glfw,btn,action,mods)

def mouseMoveCB(win_glfw,x,y):
  wmngr_glfw.motion(win_glfw,x,y)

def keyFunCB(win_glfw,key,scancode,action,mods):
  wmngr_glfw.keyinput(win_glfw,key,scancode,action,mods)

def render():
  glColor3d(1, 0, 0)
  glEnable(GL_LIGHTING)
  msh.draw()

def main():
  global msh, wmngr_glfw
  msh = dfm2.Mesh("../test_inputs/bunny_2k.ply");
  msh.scale_xyz(0.02)

  wmngr_glfw = dfm2.WindowManagerGLFW(1.0)

  glfw.init()
  win_glfw = glfw.create_window(640, 480, 'Hello World', None, None)
  glfw.make_context_current(win_glfw)

  dfm2.setSomeLighting()
  glEnable(GL_DEPTH_TEST)

  glfw.set_mouse_button_callback(win_glfw, mouseButtonCB)
  glfw.set_cursor_pos_callback(win_glfw,  mouseMoveCB)
  glfw.set_key_callback(win_glfw, keyFunCB)

  while not glfw.window_should_close(win_glfw):
    glClearColor(1, 1, 1, 1)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    wmngr_glfw.camera.set_gl_camera()
    render()
    glfw.swap_buffers(win_glfw)
    glfw.poll_events()
    if wmngr_glfw.isClose:
      break
  glfw.destroy_window(win_glfw)
  glfw.terminate()
  print("closed")

if __name__ == "__main__":
  main()