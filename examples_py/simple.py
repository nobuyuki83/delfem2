from OpenGL.GL import *
from OpenGL.GLUT import *

import numpy as np

import sys
sys.path.append("../python")
import dfm2

win = None

def display():
  global win
  glClearColor(0.3, 0.5, 0.8, 1.0)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  glEnable(GL_DEPTH_TEST)
  
  win.camera.set_gl_camera()

  glEnable(GL_LIGHTING)
  glutSolidTeapot(0.1)

  glDisable(GL_LIGHTING)
  glLineWidth(3)  
  dfm2.gl.draw_axis(size=0.2)
  glutSwapBuffers()

def reshape(width, height):
  glViewport(0, 0, width, height)


def idle():
  glutPostRedisplay()


def keyboard(bkey, x, y):
  global win
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
  global win

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

  win = dfm2.glut.WindowManager(0.3)

  dfm2.dfm2.set_some_lighting()

  ####
  # Turn the flow of control over to GLUT
  glutMainLoop()


if __name__ == "__main__":
  main()
