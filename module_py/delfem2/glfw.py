from OpenGL.GL import *
import numpy

import glfw
from .gl import *
from .libdelfem2 import *

class WindowManagerGLFW:
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
    if self.button == glfw.MOUSE_BUTTON_LEFT:
      if self.modifier == glfw.MOD_ALT:  # shift
        self.mouse_x, self.mouse_y = self.camera.rotation(x, y, self.mouse_x, self.mouse_y, win_w, win_h)
      if self.modifier == glfw.MOD_SHIFT:
        self.mouse_x, self.mouse_y = self.camera.translation(x, y, self.mouse_x, self.mouse_y, win_w, win_h)
    (x, y) = glfw.get_cursor_pos(win_glfw)
    self.mouse_x = (2.0 * x - win_w) / win_w
    self.mouse_y = (win_h - 2.0 * y) / win_h



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
    self.list_func_mouse = []
    self.list_func_motion = []
    self.list_func_step_time = []
    self.list_func_draw = []
    self.color_bg = (1,1,1)
    glEnable(GL_DEPTH_TEST)

  def draw_loop(self):
    """
    Enter the draw loop

    render -- a function to render
    """
    glfw.set_mouse_button_callback(self.win, self.mouse)
    glfw.set_cursor_pos_callback(self.win, self.motion)
    glfw.set_key_callback(self.win, self.keyinput)
    glfw.set_window_size_callback(self.win, self.window_size)
    while not glfw.window_should_close(self.win):
      glClearColor(self.color_bg[0], self.color_bg[1], self.color_bg[2], 1.0)
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
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
    mMV = glGetFloatv(GL_MODELVIEW_MATRIX)
    mPj = glGetFloatv(GL_PROJECTION_MATRIX)
    src = screenUnProjection(numpy.array([float(self.wm.mouse_x),float(self.wm.mouse_y),0.0]),
                             mMV, mPj)
    dir = screenUnProjectionDirection((0,0,1), mMV,mPj)
    for func_mouse in self.list_func_mouse:
      func_mouse(btn,action,mods, src, dir, self.wm.camera.view_height)

  def motion(self,win0,x,y):
    if self.wm.button == -1:
      return
    self.wm.motion(win0,x,y)
    mMV = glGetFloatv(GL_MODELVIEW_MATRIX)
    mPj = glGetFloatv(GL_PROJECTION_MATRIX)
    src0 = screenUnProjection(numpy.array([float(self.wm.mouse_pre_x),float(self.wm.mouse_pre_y),0.0]),
                             mMV, mPj)
    src1 = screenUnProjection(numpy.array([float(self.wm.mouse_x),float(self.wm.mouse_y),0.0]),
                             mMV, mPj)
    dir = screenUnProjectionDirection((0,0,1), mMV,mPj)
#    print(src0,src1,dir)
    for func_motion in self.list_func_motion:
      func_motion(src0,src1,dir)

  def keyinput(self,win0,key,scancode,action,mods):
    self.wm.keyinput(win0,key,scancode,action,mods)

  def window_size(self,win0,w,h):
    glViewport(0,0,w,h)



def winDraw3d(list_obj:list,
              winsize=(400,300),
              bgcolor=(1,1,1),
              glsl_vrt="",
              glsl_frg="",
              camera_eye_up=(+0.0,+0.0,-1.0, +0.0,+1.0,+0.0),
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
  if len(camera_eye_up) == 6:
    window.wm.camera.set_rotation(camera_eye_up)
  #### initalizing opengl
  setSomeLighting()
  glEnable(GL_POLYGON_OFFSET_FILL )
  glPolygonOffset( 1.1, 4.0 )
  glUseProgram(id_shader_program)
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

class DepthColorBuffer:
  def __init__(self,
    win_size:list,
    format_color:str,
    is_depth:bool):
    self.win = WindowGLFW(isVisible=False)
    self.fbm = FrameBufferManager()
    glew_init()
    self.fbm.set_buffer_size(win_size[0],win_size[1], format_color,is_depth)
  def start(self):
    self.fbm.start()
  def end(self):
    self.fbm.end()
  def close(self):
    self.end()
    self.win.close()

#def take_depth_shot(render_func, sampler:GPUSampler, buffer:DepthColorBuffer):

