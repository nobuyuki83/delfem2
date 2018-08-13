from OpenGL.GL import *
from OpenGL.GLUT import *

import numpy as np
import utility_gl


def glut_print(x, y, font, text, color):
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glMatrixMode(GL_MODELVIEW)
  glLoadIdentity()
  glColor3f(color[0], color[1], color[2])
  glRasterPos2f(x, y)
  for ch in text:
    glutBitmapCharacter(font, ctypes.c_int(ord(ch)))


def screenProjection(v0, scale, Rot, trans):
  v1 = utility_gl.mult_mat_vec_3(Rot, v0)
  v2 = [trans[0] + scale * v1[0], trans[1] - scale * v1[1], v1[2] * scale]
  return v2

def draw_sphere(pos, rad, color):
  if pos is None:
    return
  glColor3d(color[0], color[1], color[2])
  glTranslatef(+pos[0], +pos[1], +pos[2])
  glutSolidSphere(rad, 32, 32)
  glTranslatef(-pos[0], -pos[1], -pos[2])

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
  pos0 = utility_gl.rot_trans(pos0,Rot,trans)
  pos1 = utility_gl.rot_trans(pos1,Rot,trans)
  pos2 = utility_gl.rot_trans(pos2,Rot,trans)
  pos3 = utility_gl.rot_trans(pos3,Rot,trans)
  pos4 = utility_gl.rot_trans(pos4,Rot,trans)
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


def draw_rect(w0, h0, cx, cy, off_z):
  glVertex3d(-0.5 * w0 + cx, -0.5 * h0 + cy, off_z)
  glVertex3d(+0.5 * w0 + cx, -0.5 * h0 + cy, off_z)
  glVertex3d(+0.5 * w0 + cx, +0.5 * h0 + cy, off_z)
  glVertex3d(-0.5 * w0 + cx, +0.5 * h0 + cy, off_z)

camera = None
modifier = 0
mouse_x = 0.0
mouse_y = 0.0
####

def display():
  global camera
  glClearColor(0.3, 0.5, 0.8, 1.0)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  glEnable(GL_DEPTH_TEST)
  ####
  camera.set_gl_camera(600, 600)
  ####
  w0 = 0.27
  h0 = 0.19
  off_z = -0.001
  glBegin(GL_QUADS)
  glColor3d(1, 1, 1)
  draw_rect(w0, h0, 0, 0, off_z)
  #glColor3d(0, 0, 0)
  #draw_rect(0.03, 0.03, -0.12, +0.08, 0.0)
  #draw_rect(0.03, 0.03, +0.12, +0.08, 0.0)
  #draw_rect(0.03, 0.03, -0.12, -0.08, 0.0)
  #draw_rect(0.03, 0.03, +0.12, -0.08, 0.0)
  glEnd()

  glLineWidth(3)
  draw_axis(size=0.1)

  glutSwapBuffers()


def reshape(width, height):
  glViewport(0, 0, width, height)


def idle():
  glutPostRedisplay()


def keyboard(bkey, x, y):
  global edit_mode, view_mode, view_face_mode, path_csv_img
  global scale, quat, trans
  key = bkey.decode("utf-8")
  if key == 'q':
    exit()
  glutPostRedisplay()

def special(key, x, y):
  global camera
  if key == int(GLUT_KEY_PAGE_UP):
    camera.scale *= 1.03
  elif key == int(GLUT_KEY_PAGE_DOWN):
    camera.scale *= (1.0 / 1.03)
  glutPostRedisplay()


def mouse(button, state, x, y):
  global modifier, mouse_x, mouse_y
  modifier = glutGetModifiers()
  mouse_x, mouse_y = utility_gl.mouse(x, y)
  glutPostRedisplay()


def motion(x, y):
  global camera, mouse_x, mouse_y
  if modifier == 1:
    mouse_x, mouse_y, camera.quat, camera.trans = utility_gl.motion("ROT", x, y, mouse_x, mouse_y, camera.quat,
                                                                    camera.trans, camera.view_height)
  if modifier == 2:
    mouse_x, mouse_y, camera.quat, camera.trans = utility_gl.motion("TRANS", x, y, mouse_x, mouse_y, camera.quat,
                                                                    camera.trans, camera.view_height)
  glutPostRedisplay()


def main():
  global camera

  ###################33
  # GLUT Window Initialization
  glutInit()
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)  # zBuffer
  glutInitWindowSize(600, 600)
  glutInitWindowPosition(100, 100)
  glutCreateWindow("Visualizzatore_2.0")

  # Register callbacks
  glutReshapeFunc(reshape)
  glutDisplayFunc(display)
  glutMouseFunc(mouse)
  glutMotionFunc(motion)
  glutKeyboardFunc(keyboard)
  glutSpecialFunc(special)
  glutIdleFunc(idle)

  camera = utility_gl.Camera(0.3)

  ####
  # Turn the flow of control over to GLUT
  glutMainLoop()


if __name__ == "__main__":
  main()
