from OpenGL.GL import *
from OpenGL.GLUT import *

import numpy as np

import sys
sys.path.append("../python")
import dfm2

def glut_print(x, y, font, text, color):
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glMatrixMode(GL_MODELVIEW)
  glLoadIdentity()
  glColor3f(color[0], color[1], color[2])
  glRasterPos2f(x, y)
  for ch in text:
    glutBitmapCharacter(font, ctypes.c_int(ord(ch)))

def draw_sphere(pos, rad, color):
  if pos is None:
    return
  glColor3d(color[0], color[1], color[2])
  glTranslatef(+pos[0], +pos[1], +pos[2])
  glutSolidSphere(rad, 32, 32)
  glTranslatef(-pos[0], -pos[1], -pos[2])

win = None
mshtri = None
####

def display():
  global camera
  glClearColor(0.3, 0.5, 0.8, 1.0)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  glEnable(GL_DEPTH_TEST)
  ####
  win.camera.set_gl_camera(600, 600)
  ####
  glColor3d(1, 0, 0)

  glEnable(GL_LIGHTING)
  mshtri.draw()

  glDisable(GL_LIGHTING)
  glLineWidth(3)
  dfm2.gl.draw_axis(size=0.5)
  glutSwapBuffers()


def reshape(width, height):
  glViewport(0, 0, width, height)


def idle():
  glutPostRedisplay()


def keyboard(bkey, x, y):
  key = bkey.decode("utf-8")
  if key == 'q':
    exit()
  glutPostRedisplay()

def special(key, x, y):
  global win
  win.special(key,x,y)

def mouse(button, state, x, y):
  global win
  win.mouse(button,state,x,y)

def motion(x, y):
  global win
  win.motion(x,y)


def main():
  global win, mshtri

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

  win = dfm2.glut.WindowManager(1.0)
  dfm2.dfm2.set_some_lighitng()

  mshtri = dfm2.dfm2.MeshTri();
  mshtri.read("../test_inputs/bunny_2k.ply")
  mshtri.scale_xyz(0.01)

  ####
  # Turn the flow of control over to GLUT
  glutMainLoop()


if __name__ == "__main__":
  main()
