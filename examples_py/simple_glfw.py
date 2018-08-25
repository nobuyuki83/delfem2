from OpenGL.GL import *
import glfw

import sys
sys.path.append("../python")
import dfm2

mshelm = None
isClose = False

camera = None
modifier = 0
mouse_x = 0
mouse_y = 0

def mouseButtonCB(win_glfw,button,action,mods):
  global mouse_x,mouse_y,modifier
  (win_w,win_h) = glfw.get_window_size(win_glfw)
  (x,y) = glfw.get_cursor_pos(win_glfw)
  mouse_x = (2.0 * x - win_w) / win_w
  mouse_y = (win_h - 2.0 * y) / win_h
  modifier = mods
  print(win_w, win_h, (x, y), (mouse_x,mouse_y), mods)

def mouseMoveCB(win_glfw,x,y):
  global mouse_x,mouse_y
  (win_w,win_h) = glfw.get_window_size(win_glfw)
  if modifier == glfw.MOD_ALT:  # shift
    mouse_x, mouse_y = camera.rotation(x,y,mouse_x,mouse_y,win_w,win_h)
  if modifier == glfw.MOD_SHIFT:
    mouse_x, mouse_y = camera.translation(x, y, mouse_x, mouse_y, win_w, win_h)
  print(mouse_x,mouse_y,x,y)

def keyFunCB(win_glfw,key,scancode,action,mods):
  global isClose,camera
  print(key)
  if key == glfw.KEY_Q and action == glfw.PRESS:
    isClose = True
  if key == glfw.KEY_PAGE_UP:
    camera.scale *= 1.03
  if key == glfw.KEY_PAGE_DOWN:
    camera.scale /= 1.03

def render():
  glColor3d(1, 0, 0)
  glEnable(GL_LIGHTING)
  mshelm.drawFace_elemWiseNorm()

def main():
  global mshelm,camera
  mshelm = dfm2.dfm2.MeshElem("../test_inputs/bunny_2k.ply");
  mshelm.scaleXYZ(0.02)

  camera = dfm2.gl.Camera(1.0)

  glfw.init()
  win_glfw = glfw.create_window(640, 480, 'Hello World', None, None)
  glfw.make_context_current(win_glfw)

  dfm2.setSomeLighting()
  glEnable(GL_DEPTH_TEST)

#  glfwSetWindowSizeCallback(window, windowSizeCB);
  glfw.set_mouse_button_callback(win_glfw, mouseButtonCB)
  glfw.set_cursor_pos_callback(win_glfw,  mouseMoveCB)
  glfw.set_key_callback(win_glfw, keyFunCB)

  #glfw.set_char_callback(win_glfw, charFunCB);
#  glfwSetDropCallback(window, dropCB);
  #  glfwSetScrollCallback(window, mouseScrollCB);

  while not glfw.window_should_close(win_glfw):
    glClearColor(1, 1, 1, 1)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    camera.set_gl_camera()
    render()
    glfw.swap_buffers(win_glfw)
    glfw.poll_events()
    if isClose:
      break
  glfw.destroy_window(win_glfw)
  glfw.terminate()

  print("closed")

if __name__ == "__main__":
  main()