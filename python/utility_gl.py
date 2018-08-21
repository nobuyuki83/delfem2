import math
from OpenGL.GL import *


def draw_axis(size=1.0):
  glBegin(GL_LINES)
  glColor3d(1,0,0)
  glVertex3d(0,0,0)
  glVertex3d(size,0,0)
  ####
  glColor3d(0,1,0)
  glVertex3d(0,0,0)
  glVertex3d(0,size,0)
  ####
  glColor3d(0,0,1)
  glVertex3d(0,0,0)
  glVertex3d(0,0,size)
  glEnd()

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


def rot_matrix_cartesian(vec):
  sqt = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]
  mat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
  if sqt < 1.0e-20:  # infinitesmal rotation approximation
    mat[0][0] = 1
    mat[0][1] = -vec[2]
    mat[0][2] = +vec[1]
    mat[1][0] = +vec[2]
    mat[1][1] = 1
    mat[1][2] = -vec[0]
    mat[2][0] = -vec[1]
    mat[2][1] = +vec[0]
    mat[2][2] = 1
    return mat
  t = math.sqrt(sqt)
  invt = 1.0 / t
  n = [vec[0] * invt, vec[1] * invt, vec[2] * invt]
  c0 = math.cos(t)
  s0 = math.sin(t)
  mat[0][0] = c0 + (1 - c0) * n[0] * n[0]
  mat[0][1] = -n[2] * s0 + (1 - c0) * n[0] * n[1]
  mat[0][2] = +n[1] * s0 + (1 - c0) * n[0] * n[2]
  mat[1][0] = +n[2] * s0 + (1 - c0) * n[1] * n[0]
  mat[1][1] = c0 + (1 - c0) * n[1] * n[1]
  mat[1][2] = -n[0] * s0 + (1 - c0) * n[1] * n[2]
  mat[2][0] = -n[1] * s0 + (1 - c0) * n[2] * n[0]
  mat[2][1] = +n[0] * s0 + (1 - c0) * n[2] * n[1]
  mat[2][2] = c0 + (1 - c0) * n[2] * n[2]
  return mat

def add_scaled_3(vec0,vec1,s):
  vec2 = [0,0,0]
  vec2[0] = vec0[0] + s*vec1[0]
  vec2[1] = vec0[1] + s*vec1[1]
  vec2[2] = vec0[2] + s*vec1[2]
  return vec2


def get_quaternion_rot_matrix(mat):
  smat = [
    1 + mat[0][0] + mat[1][1] + mat[2][2],
    mat[2][1] - mat[1][2],
    mat[0][2] - mat[2][0],
    mat[1][0] - mat[0][1],
    mat[2][1] - mat[1][2],
    1 + mat[0][0] - mat[1][1] - mat[2][2],
    mat[0][1] + mat[1][0],
    mat[0][2] + mat[2][0],
    mat[0][2] - mat[2][0],
    mat[1][0] + mat[0][1],
    1 - mat[0][0] + mat[1][1] - mat[2][2],
    mat[1][2] + mat[2][1],
    mat[1][0] - mat[0][1],
    mat[0][2] + mat[2][0],
    mat[1][2] + mat[2][1],
    1 - mat[0][0] - mat[1][1] + mat[2][2]]

  imax = 0
  if smat[1 * 4 + 1] > smat[imax * 4 + imax]:
    imax = 1
  if smat[2 * 4 + 2] > smat[imax * 4 + imax]:
    imax = 2
  if smat[3 * 4 + 3] > smat[imax * 4 + imax]:
    imax = 3

  eparam2 = [0, 0, 0, 0]
  eparam2[imax] = 0.5 * math.sqrt(smat[imax * 4 + imax])
  for k in range(4):
    if k == imax:
      continue
    eparam2[k] = smat[imax * 4 + k] * 0.25 / eparam2[imax]
  return eparam2


def rot_matrix_quaternion(q):
  x2 = q[1] * q[1] * 2.0
  y2 = q[2] * q[2] * 2.0
  z2 = q[3] * q[3] * 2.0
  xy = q[1] * q[2] * 2.0
  yz = q[2] * q[3] * 2.0
  zx = q[3] * q[1] * 2.0
  xw = q[1] * q[0] * 2.0
  yw = q[2] * q[0] * 2.0
  zw = q[3] * q[0] * 2.0
  R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
  R[0][0] = 1.0 - y2 - z2
  R[1][0] = xy + zw
  R[2][0] = zx - yw
  R[0][1] = xy - zw
  R[1][1] = 1.0 - z2 - x2
  R[2][1] = yz + xw
  R[0][2] = zx + yw
  R[1][2] = yz - xw
  R[2][2] = 1.0 - x2 - y2
  return R


def mult_GlAffineMatrix(m, p):
  v = [0, 0, 0]
  v[0] = m[0][0] * p[0] + m[1][0] * p[1] + m[2][0] * p[2] + m[3][0]
  v[1] = m[0][1] * p[0] + m[1][1] * p[1] + m[2][1] * p[2] + m[3][1]
  v[2] = m[0][2] * p[0] + m[1][2] * p[1] + m[2][2] * p[2] + m[3][2]
  return v


def mult_mat_vec_3(m, v):
  r = [0, 0, 0]
  r[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2]
  r[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2]
  r[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
  return r


def mat_trans_3(m):
  n = [[m[0][0],m[1][0],m[2][0]],
       [m[0][1],m[1][1],m[2][1]],
       [m[0][2],m[1][2],m[2][2]]]
  return n


def mult_mat3_mat3(a, b):
  c = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
  for i in range(3):
    for j in range(3):
      c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]
  return c

def distance_3(a,b):
  s0 = (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2])
  return math.sqrt(s0)


def rot_trans(v,R,t):
  v0 = v
  if R is not None:
    v0 = mult_mat_vec_3(R,v)
  v1 = v0
  if t is not None:
    v1 = add_scaled_3(v0,t,+1.0)
  return v1


def scale_vec3(a,d):
  v = [a[0]*d,a[1]*d,a[2]*d]
  return v


def average_vec3(list_pos):
  ave = [0,0,0]
  icnt = 0
  for pos in list_pos:
    if pos is None:
      continue
    ave[0] += pos[0]
    ave[1] += pos[1]
    ave[2] += pos[2]
    icnt += 1
  if icnt == 0:
    return None
  else:
    ave[0] /= icnt
    ave[1] /= icnt
    ave[2] /= icnt
  return ave

def std_vec3(list_pos,ave):
  sum_dist = 0
  icnt = 0
  for pos in list_pos:
    if pos is None:
      continue
    dist = distance_3(pos,ave)
    sum_dist += dist*dist
    icnt += 1
  if icnt == 0:
    return None
  sum_dist /= icnt
  return math.sqrt(sum_dist)

def affine_matrix_quaternion(q):
  x2 = q[1] * q[1] * 2.0
  y2 = q[2] * q[2] * 2.0
  z2 = q[3] * q[3] * 2.0
  xy = q[1] * q[2] * 2.0
  yz = q[2] * q[3] * 2.0
  zx = q[3] * q[1] * 2.0
  xw = q[1] * q[0] * 2.0
  yw = q[2] * q[0] * 2.0
  zw = q[3] * q[0] * 2.0

  r = [None] * 16  # empty list
  r[0] = 1.0 - y2 - z2
  r[1] = xy + zw
  r[2] = zx - yw
  r[4] = xy - zw
  r[5] = 1.0 - z2 - x2
  r[6] = yz + xw
  r[8] = zx + yw
  r[9] = yz - xw
  r[10] = 1.0 - x2 - y2
  r[3] = r[7] = r[11] = r[12] = r[13] = r[14] = 0.0
  r[15] = 1.0
  return r


def QuatMult(p, q):
  r = [None] * 4
  r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3]
  r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2]
  r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1]
  r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0]
  return r


def motion(edit_mode, x, y, mouse_x, mouse_y, quat, trans, view_height):
  assert len(trans) == 2
  assert len(quat) == 4
  viewport = glGetIntegerv(GL_VIEWPORT)
  win_w = viewport[2]
  win_h = viewport[3]
  mov_end_x = (2.0 * x - win_w) / win_w
  mov_end_y = (win_h - 2.0 * y) / win_h
  dx = mov_end_x - mouse_x
  dy = mov_end_y - mouse_y
  if edit_mode == "ROT":
    a = math.sqrt(dx * dx + dy * dy)
    ar = a * 0.5
    dq = [math.cos(ar), -dy * math.sin(ar) / a, dx * math.sin(ar) / a, 0.0]
    if a != 0.0:
      quat = QuatMult(dq, quat)
  elif edit_mode == "TRANS":
    trans[0] += dx * view_height * 0.5
    trans[1] += dy * view_height * 0.5

  mouse_x = mov_end_x
  mouse_y = mov_end_y
  return mouse_x, mouse_y, quat, trans


def mouse(x, y):
  viewport = glGetIntegerv(GL_VIEWPORT)
  win_w = viewport[2]
  win_h = viewport[3]
  mouse_x = (2.0 * x - win_w) / win_w
  mouse_y = (win_h - 2.0 * y) / win_h
  return mouse_x, mouse_y


class Camera:
  def __init__(self, view_height):
    self.view_height = view_height
    self.scale = 1.0
    self.trans = [0, 0]
    self.quat = [1, 0, 0, 0]
    self.fovy = 60  # degree

  def set_gl_camera(self, win_w, win_h):
    depth = self.view_height / (self.scale * math.tan(0.5 * self.fovy * 3.1415 / 180.0))
    asp = float(win_w) / win_h
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(-self.view_height / self.scale * asp,
            +self.view_height / self.scale * asp,
            -self.view_height / self.scale,
            +self.view_height / self.scale,
            -depth * 10,
            +depth * 10)

    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    glTranslated(self.trans[0], self.trans[1], -depth)

    Rview = affine_matrix_quaternion(self.quat)
    glMultMatrixd(Rview)
