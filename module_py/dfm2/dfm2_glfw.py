from OpenGL.GL import *
import numpy

import glfw
from .dfm2_gl import Camera
from .dfm2 import *

class WindowManagerGLFW:
  def __init__(self,view_height):
    self.camera = Camera(view_height)
    self.modifier = 0
    self.mouse_x = 0.0
    self.mouse_y = 0.0
    self.button = -1
    self.isClose = False

  def keyinput(self,win_glfw,key,scancode,action,mods):
    if key == glfw.KEY_Q and action == glfw.PRESS:
      self.isClose = True
    if key == glfw.KEY_PAGE_UP:
      self.camera.scale *= 1.03
    if key == glfw.KEY_PAGE_DOWN:
      self.camera.scale /= 1.03

  def mouse(self,win_glfw,btn,action,mods):
    (win_w, win_h) = glfw.get_window_size(win_glfw)
    (x, y) = glfw.get_cursor_pos(win_glfw)
    self.mouse_x = (2.0 * x - win_w) / win_w
    self.mouse_y = (win_h - 2.0 * y) / win_h
    self.modifier = mods
    if action == glfw.PRESS:
      self.button = btn
    elif action == glfw.RELEASE:
      self.button = -1

  def motion(self,win_glfw, x, y):
    (win_w, win_h) = glfw.get_window_size(win_glfw)
    if self.button == glfw.MOUSE_BUTTON_LEFT:
      if self.modifier == glfw.MOD_ALT:  # shift
        self.mouse_x, self.mouse_y = self.camera.rotation(x, y, self.mouse_x, self.mouse_y, win_w, win_h)
      if self.modifier == glfw.MOD_SHIFT:
        self.mouse_x, self.mouse_y = self.camera.translation(x, y, self.mouse_x, self.mouse_y, win_w, win_h)


class WindowGLFW:
  """
  class to manage the glfw window
  """
  def __init__(self,view_height=1.0,winsize=(400,300),isVisible=True):
    if glfw.init() == GL_FALSE:
      print("GLFW couldn't not initialize!")
    if not isVisible:
      glfw.window_hint(glfw.VISIBLE, False)
    self.win = glfw.create_window(winsize[0], winsize[1], '3D Window', None, None)
    glfw.make_context_current(self.win)
    ###
    self.wm = WindowManagerGLFW(view_height)
    self.draw_func = None
    glEnable(GL_DEPTH_TEST)

  def draw_loop(self,list_draw_func):
    """
    Enter the draw loop

    render -- a function to render
    """
    glfw.set_mouse_button_callback(self.win, self.mouse)
    glfw.set_cursor_pos_callback(self.win, self.motion)
    glfw.set_key_callback(self.win, self.keyinput)
    while not glfw.window_should_close(self.win):
      glClearColor(1, 1, 1, 1)
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
      self.wm.camera.set_gl_camera()
      for draw_func in list_draw_func:
        draw_func()
      glfw.swap_buffers(self.win)
      glfw.poll_events()
      if self.wm.isClose:
        break
    glfw.destroy_window(self.win)
    glfw.terminate()

  def mouse(self, win0,btn,action,mods):
    self.wm.mouse(win0,btn,action,mods)

  def motion(self,win0,x,y):
    self.wm.motion(win0,x,y)

  def keyinput(self,win0,key,scancode,action,mods):
    self.wm.keyinput(win0,key,scancode,action,mods)



def winDraw3d(list_obj,winsize=(400,300)):
  """
  draw the input object into openGL window

  obj -- the object to draw.
    this object need to have a method "draw()"

  winsize -- the size of the window (width,height)
  """
  #### initialize window
  window = WindowGLFW(winsize=winsize)
  #### adjust scale
  aabb3 = AABB3()
  for obj in list_obj:
    aabb3.add_minmax_xyz(obj.minmax_xyz())
  if not aabb3.isActive:
    aabb3.set_minmax_xyz(-1,+1, -1,+1, -1,+1)
  window.wm.camera.adjust_scale_trans(aabb3.list_xyz())
  #### initalizing opengl
  setSomeLighting()  
  glEnable(GL_POLYGON_OFFSET_FILL )
  glPolygonOffset( 1.1, 4.0 )
  #### enter loop
  window.draw_loop([x.draw for x in list_obj])


def imgDraw3d(list_obj,winsize=(400,300)):
  """
  draw the input object into Numpy uint8 array

  obj -- the object to draw

  winsize -- the size of the window
  """
  #### initialize window
  window = WindowGLFW(1.0,winsize=winsize,isVisible=False)
  #### set camera
  aabb3 = AABB3()
  for obj in list_obj:
    aabb3.add_minmax(obj.minmax_xyz())
  window.wm.camera.adjust_scale_trans(aabb3.list_xyz())
  #### initialize opengl
  setSomeLighting()
  glEnable(GL_POLYGON_OFFSET_FILL )
  glPolygonOffset( 1.1, 4.0 )
  #### draw
  glClearColor(1, 1, 1, 1)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  window.wm.camera.set_gl_camera()
  for obj in list_obj:
    obj.draw()
  glFlush()
  #### read pixel to img
  glPixelStorei(GL_PACK_ALIGNMENT, 1)
  (buff_w,buff_h) = glGetIntegerv(GL_VIEWPORT)[2:]
  bytes_img = glReadPixels(0, 0, buff_w, buff_h, GL_RGB, GL_UNSIGNED_BYTE)
  img = numpy.frombuffer(bytes_img, dtype=numpy.uint8)
  #### terminate window
  glfw.destroy_window(window.win)
  glfw.terminate()
  #### reshape the img array
  img = numpy.reshape(img,(buff_h,buff_w,3))
  return img