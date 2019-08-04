####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import OpenGL.GL as gl
import glfw

import numpy
from typing import List

from ._gl import Camera, screenUnProjection, screenUnProjectionDirection, AxisXYZ

from .c_gl import CppFrameBufferManager
from .c_gl import glew_init, setSomeLighting

from ..c_core import AABB3


class NavigationGLFW:
  def __init__(self,view_height):
    self.camera = Camera(view_height)
    self.modifier = 0
    self.mouse_x = 0.0
    self.mouse_y = 0.0
    self.mouse_pre_x = 0.0
    self.mouse_pre_y = 0.0
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
    self.mouse_pre_x,self.mouse_pre_y = self.mouse_x, self.mouse_y
    (x, y) = glfw.get_cursor_pos(win_glfw)
    self.mouse_x = (2.0 * x - win_w) / win_w
    self.mouse_y = (win_h - 2.0 * y) / win_h
    if self.button == glfw.MOUSE_BUTTON_LEFT:
      if self.modifier == glfw.MOD_ALT:  # shift
        self.camera.rotation(self.mouse_x, self.mouse_y, self.mouse_pre_x, self.mouse_pre_y)
      if self.modifier == glfw.MOD_SHIFT:
        self.camera.translation(self.mouse_x, self.mouse_y, self.mouse_pre_x, self.mouse_pre_y)



class WindowGLFW:
  """
  class to manage the glfw window
  """
  def __init__(self,view_height=1.0,winsize=(400,300),isVisible=True):
    if glfw.init() == gl.GL_FALSE:
      print("GLFW couldn't not initialize!")
    if not isVisible:
      glfw.window_hint(glfw.VISIBLE, False)
    self.win = glfw.create_window(winsize[0], winsize[1], '3D Window', None, None)
    glfw.make_context_current(self.win)
    ###
    self.wm = NavigationGLFW(view_height)
    self.list_func_mouse = []
    self.list_func_motion = []
    self.list_func_step_time = []
    self.list_func_draw = []
    self.color_bg = (1,1,1)
    gl.glEnable(gl.GL_DEPTH_TEST)

  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc_val, exc_tb):
#    self.close()
    pass

  def draw_loop(self):
    """
    Enter the draw loop

    render -- a function to render
    """
    glfw.set_mouse_button_callback(self.win, self.mouse)
    glfw.set_cursor_pos_callback(self.win, self.motion)
    glfw.set_key_callback(self.win, self.keyinput)
#    glfw.set_window_size_callback(self.win, self.window_size)
    while not glfw.window_should_close(self.win):
      gl.glClearColor(self.color_bg[0], self.color_bg[1], self.color_bg[2], 1.0)
      gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
      gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
      gl.glPolygonOffset(1.1, 4.0)
      self.wm.camera.set_gl_camera()
      for func_step_time in self.list_func_step_time:
        func_step_time()
      for draw_func in self.list_func_draw:
        draw_func()
      glfw.swap_buffers(self.win)
      glfw.poll_events()
      if self.wm.isClose:
        break
    self.close()

  def close(self):
    glfw.destroy_window(self.win)
    glfw.terminate()

  def mouse(self, win0,btn,action,mods):
    self.wm.mouse(win0,btn,action,mods)
    mMV = gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX)
    mPj = gl.glGetFloatv(gl.GL_PROJECTION_MATRIX)
    src = screenUnProjection(numpy.array([float(self.wm.mouse_x),float(self.wm.mouse_y),0.0]),
                             mMV, mPj)
    dir = screenUnProjectionDirection(numpy.array([0,0,1]), mMV,mPj)
    for func_mouse in self.list_func_mouse:
      func_mouse(btn,action,mods, src, dir, self.wm.camera.view_height)

  def motion(self,win0,x,y):
    if self.wm.button == -1:
      return
    self.wm.motion(win0,x,y)
    mMV = gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX)
    mPj = gl.glGetFloatv(gl.GL_PROJECTION_MATRIX)
    src0 = screenUnProjection(numpy.array([float(self.wm.mouse_pre_x),float(self.wm.mouse_pre_y),0.0]),
                             mMV, mPj)
    src1 = screenUnProjection(numpy.array([float(self.wm.mouse_x),float(self.wm.mouse_y),0.0]),
                             mMV, mPj)
    dir = screenUnProjectionDirection(numpy.array([0,0,1]), mMV,mPj)
    for func_motion in self.list_func_motion:
      func_motion(src0,src1,dir)

  def keyinput(self,win0,key,scancode,action,mods):
    self.wm.keyinput(win0,key,scancode,action,mods)

#  def window_size(self,win0,w,h):
#    gl.glViewport(0,0,w,h) # because of multisampling this doesn't work



def winDraw3d(list_obj:list,
              winsize=(400,300),
              bgcolor=(1,1,1),
              glsl_vrt="",
              glsl_frg="",
              camera_orientation=(+0.0,+0.0,-1.0, +0.0,+1.0,+0.0),
              camera_scale=1.0):
  """
  draw the input object into openGL window

  obj -- the object to draw.
    this object need to have a method "draw()"

  winsize -- the size of the window (width,height)
  """
  #### initialize window
  window = WindowGLFW(winsize=winsize)
  window.color_bg = bgcolor
  for obj in list_obj:
    if hasattr(obj, 'init_gl'):
      obj.init_gl()
    if hasattr(obj, 'mouse'):
      window.list_func_mouse.append(obj.mouse)
    if hasattr(obj, 'motion'):
      window.list_func_motion.append(obj.motion)
    if hasattr(obj, "draw"):
      window.list_func_draw.append(obj.draw)
    if hasattr(obj, "step_time"):
      window.list_func_step_time.append(obj.step_time)
  #### glsl compile
  id_shader_program = 0
  if glsl_vrt != "" and glsl_frg != "":
    glew_init()
    id_shader_program = setup_glsl(glsl_vrt,glsl_frg)
  #### adjust scale
  aabb3 = AABB3()
  for obj in list_obj:
    if hasattr(obj, 'minmax_xyz'):
      aabb3.add_minmax_xyz(obj.minmax_xyz())
  if not aabb3.isActive:
    aabb3.set_minmax_xyz(-1,+1, -1,+1, -1,+1)
  window.wm.camera.adjust_scale_trans(aabb3.list_xyz())
  window.wm.camera.scale = camera_scale
  #### set camera rotation
  if len(camera_orientation) == 6:
    window.wm.camera.set_rotation(camera_orientation)
  #### initalizing opengl
  setSomeLighting()
  gl.glEnable(gl.GL_POLYGON_OFFSET_FILL )
  gl.glPolygonOffset( 1.1, 4.0 )
  gl.glUseProgram(id_shader_program)
  #### enter loop
  window.draw_loop()


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
    aabb3.add_minmax_xyz(obj.minmax_xyz())
  window.wm.camera.adjust_scale_trans(aabb3.list_xyz())
  #### initialize opengl
  setSomeLighting()
  gl.glEnable(gl.GL_POLYGON_OFFSET_FILL )
  gl.glPolygonOffset( 1.1, 4.0 )
  #### draw
  gl.glClearColor(1, 1, 1, 1)
  gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
  window.wm.camera.set_gl_camera()
  for obj in list_obj:
    obj.draw()
  gl.glFlush()
  #### read pixel to img
  gl.glPixelStorei(gl.GL_PACK_ALIGNMENT, 1)
  (buff_w,buff_h) = gl.glGetIntegerv(gl.GL_VIEWPORT)[2:]
  bytes_img = gl.glReadPixels(0, 0, buff_w, buff_h, gl.GL_RGB, gl.GL_UNSIGNED_BYTE)
  img = numpy.frombuffer(bytes_img, dtype=numpy.uint8)
  #### terminate window
  glfw.destroy_window(window.win)
  glfw.terminate()
  #### reshape the img array
  img = numpy.reshape(img,(buff_h,buff_w,3))
  return img


class GPUSamplerBufferGLFW:
  def __init__(self,
               win_size: List[int],
               format_color: str,
               is_depth: bool):
    self.win = WindowGLFW(isVisible=False)
    glew_init()
    self.fbm = CppFrameBufferManager()
    self.fbm.set_buffer_size(win_size[0],win_size[1], format_color,is_depth)
    self.fbm.start()
    self.is_open = True

  def __enter__(self):
    return self

  def __exit__(self, ex_type, ex_value, trace):
    self.close()

  def __del__(self):
    self.close()

  def close(self):
    if self.is_open:
      self.fbm.end()
      self.win.close()
    self.is_open = False


