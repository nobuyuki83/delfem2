####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import math, numpy
import OpenGL.GL as gl
import sys
sys.path.append("..")

from ..util import minMaxLoc, v3_scale, v3_normalize, v3_cross
from ..util import motion_rot, motion_trans
from ..util import get_quaternion_rot_matrix, affine_matrix_quaternion
from ..util import screenUnProjection, screenUnProjectionDirection

def getOpenglInfo():
  info = """
    Vendor: {0}
    Renderer: {1}
    OpenGL Version: {2}
    Shader Version: {3}
  """.format(
    gl.glGetString(gl.GL_VENDOR),
    gl.glGetString(gl.GL_RENDERER),
    gl.glGetString(gl.GL_VERSION),
    gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION)
  )
  return info


class AxisXYZ():
  def __init__(self, size=1.0, pos=(0.0,0.0,0.0), line_width=1):
    self.size = size
    self.pos = pos
    self.line_width=line_width

  def draw(self):
    gl.glDisable(gl.GL_LIGHTING)
    gl.glLineWidth(self.line_width)
    gl.glTranslated(+self.pos[0],+self.pos[1],+self.pos[2])
    ####
    gl.glBegin(gl.GL_LINES)
    gl.glColor3d(1,0,0)
    gl.glVertex3d(0,0,0)
    gl.glVertex3d(self.size,0,0)
    ####
    gl.glColor3d(0,1,0)
    gl.glVertex3d(0,0,0)
    gl.glVertex3d(0,self.size,0)
    ####
    gl.glColor3d(0,0,1)
    gl.glVertex3d(0,0,0)
    gl.glVertex3d(0,0,self.size)
    gl.glEnd()
    ####
    gl.glTranslated(-self.pos[0],-self.pos[1],-self.pos[2])

  def minmax_xyz(self):
    return [0,self.size, 0, self.size, 0, self.size]



def draw_pyramid(lenWh, lenHh, lenZ, Rot, trans):
  pos0 = [-lenWh,-lenHh,+lenZ]
  pos1 = [+lenWh,-lenHh,+lenZ]
  pos2 = [+lenWh,+lenHh,+lenZ]
  pos3 = [-lenWh,+lenHh,+lenZ]
  pos4 = [0.0, 0.0, 0.0]
  pos0 = rot_trans(pos0,Rot,trans)
  pos1 = rot_trans(pos1,Rot,trans)
  pos2 = rot_trans(pos2,Rot,trans)
  pos3 = rot_trans(pos3,Rot,trans)
  pos4 = rot_trans(pos4,Rot,trans)
  glBegin(GL_LINES)
  glVertex3dv(pos0)
  glVertex3dv(pos1)
  glVertex3dv(pos1)
  glVertex3dv(pos2)
  glVertex3dv(pos2)
  glVertex3dv(pos3)
  glVertex3dv(pos3)
  glVertex3dv(pos0)
  glVertex3dv(pos0)
  glVertex3dv(pos4)
  glVertex3dv(pos1)
  glVertex3dv(pos4)
  glVertex3dv(pos2)
  glVertex3dv(pos4)
  glVertex3dv(pos3)
  glVertex3dv(pos4)
  glEnd()

class Camera:
  def __init__(self, view_height=1.0):
    self.view_height = view_height
    self.scale = 1.0
    self.scr_trans = [0., 0.] # position of the pivot in the screen
    self.pivot = [0., 0., 0.] # pivot location
    self.quat = [1, 0, 0, 0]
    self.fovy = 60  # degree

  def set_gl_camera(self):
    viewport = gl.glGetIntegerv(gl.GL_VIEWPORT)
    win_w = viewport[2]
    win_h = viewport[3]
    depth = self.view_height / (self.scale * math.tan(0.5 * self.fovy * 3.1415 / 180.0))
    asp = float(win_w) / win_h
    ####
    gl.glMatrixMode(gl.GL_PROJECTION)
    gl.glLoadIdentity()
    gl.glOrtho(-self.view_height / self.scale * asp,
              +self.view_height / self.scale * asp,
              -self.view_height / self.scale,
              +self.view_height / self.scale,
              -depth * 10,
              +depth * 10)

    ####
    gl.glMatrixMode(gl.GL_MODELVIEW)
    gl.glLoadIdentity()
    gl.glTranslated(self.scr_trans[0], self.scr_trans[1], -depth)
    Rview = affine_matrix_quaternion(self.quat)
    gl.glMultMatrixd(Rview)
    gl.glTranslated(self.pivot[0], self.pivot[1], self.pivot[2])

  def set_rotation(self, camera_eye_up):
    vz = camera_eye_up[0:3]
    vy = camera_eye_up[3:6]
    vz = v3_scale(v3_normalize(vz),-1.0)
    vy = v3_normalize(vy)
    vx = v3_cross(vy,vz)
    mat = [vx,vy,vz]
    self.quat = get_quaternion_rot_matrix(mat)

  def rotation(self,sx1,sy1,sx0,sy0):
    self.quat, self.trans = motion_rot(
      sx1, sy1, sx0, sy0,
      self.quat, self.scr_trans)
    return

  def translation(self,sx1,sy1,sx0,sy0):
    self.quat, self.trans = motion_trans(
      sx1, sy1, sx0, sy0, self.quat,
      self.scr_trans, self.view_height,
      self.scale)
    return

  def adjust_scale_trans(self, aPos):
    minmax_x = minMaxLoc(aPos, [1., 0., 0.])
    minmax_y = minMaxLoc(aPos, [0., 1., 0.])
    minmax_z = minMaxLoc(aPos, [0., 0., 1.])
    (win_w,win_h) = gl.glGetIntegerv(gl.GL_VIEWPORT)[2:]
    asp = float(win_w) / win_h
    vh1 = (minmax_x[1]-minmax_x[0])/asp
    vh0 = (minmax_y[1]-minmax_y[0])
    self.pivot[0] = -0.5*(minmax_x[0]+minmax_x[1])
    self.pivot[1] = -0.5*(minmax_y[0]+minmax_y[1])
    self.pivot[2] = -0.5*(minmax_z[0]+minmax_z[1])
    self.scr_trans[0] = 0.0
    self.scr_trans[1] = 0.0
    self.view_height = max(vh0,vh1)
    self.scale = 1.0
