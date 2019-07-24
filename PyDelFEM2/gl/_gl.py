####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import math, numpy
import OpenGL.GL as gl

from ..util import minMaxLoc, v3_scale, v3_normalize, v3_cross
from ..util import motion_rot, motion_trans
from ..util import get_quaternion_rot_matrix, affine_matrix_quaternion
from ..util import screenUnProjection, screenUnProjectionDirection

from enum import Enum

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


'''
def glhLookAtf2(matrix:List[float], eyePosition3D:List[float], 
                  float *center3D, float *upVector3D )
{
   float forward[3], side[3], up[3];
   float matrix2[16], resultMatrix[16];
   // --------------------
   forward[0] = center3D[0] - eyePosition3D[0];
   forward[1] = center3D[1] - eyePosition3D[1];
   forward[2] = center3D[2] - eyePosition3D[2];
   NormalizeVector(forward);
   // --------------------
   // Side = forward x up
   ComputeNormalOfPlane(side, forward, upVector3D);
   NormalizeVector(side);
   --------------------
   // Recompute up as: up = side x forward
   ComputeNormalOfPlane(up, side, forward);
   // --------------------
   matrix2[0] = side[0];
   matrix2[4] = side[1];
   matrix2[8] = side[2];
   matrix2[12] = 0.0;
   // --------------------
   matrix2[1] = up[0];
   matrix2[5] = up[1];
   matrix2[9] = up[2];
   matrix2[13] = 0.0;
   // --------------------
'''

class CAMERA_ROT_MODE(Enum):
  YTOP = 1
  ZTOP = 2
  TBALL = 3


class Camera:
  def __init__(self, view_height=1.0):
    self.view_height = view_height
    self.scale = 1.0
    self.scr_trans = [0., 0.] # position of the pivot in the screen
    self.pivot = [0., 0., 0.] # pivot location
    self.quat = [1, 0, 0, 0]
    self.camera_rot_mode = CAMERA_ROT_MODE.TBALL
    self.fovy = 60  # degree
    self.theta = 0.0
    self.psi = 0.0

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
    if self.camera_rot_mode == CAMERA_ROT_MODE.TBALL:
      Rview = affine_matrix_quaternion(self.quat)
      gl.glMultMatrixd(Rview)
    elif self.camera_rot_mode == CAMERA_ROT_MODE.YTOP:
      x = math.sin(self.theta)*math.cos(self.psi)
      z = math.cos(self.theta)*math.cos(self.psi)
      y = math.sin(self.psi)
      ey = numpy.array([x,y,z])
      up = numpy.array([0,1,0])
      up =  up-ey*numpy.dot(ey,up)
      up /= numpy.linalg.norm(up)
      vx = v3_normalize(v3_cross(up, ey))
      mMV = [vx[0], up[0], ey[0], 0,
             vx[1], up[1], ey[1], 0,
             vx[2], up[2], ey[2], 0,
             0,     0,     0,     1]
      mMV = numpy.array(mMV)
      gl.glMultMatrixf(mMV)
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
    if self.camera_rot_mode == CAMERA_ROT_MODE.TBALL:
      self.quat, self.trans = motion_rot(
        sx1, sy1, sx0, sy0,
        self.quat, self.scr_trans)
    elif self.camera_rot_mode == CAMERA_ROT_MODE.YTOP:
      self.theta += sx1-sx0
      self.psi   += sy1-sy0
      self.psi = min(max(-math.pi*0.499,self.psi),math.pi*0.499)
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
