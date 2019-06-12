from OpenGL.GL import *
from OpenGL.GLUT import *

import sys
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glut

wmngr_glut = None

def display():
  global wmngr_glut
  glClearColor(0.3, 0.5, 0.8, 1.0)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  glEnable(GL_DEPTH_TEST)
  
  wmngr_glut.camera.set_gl_camera()

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
  global wmngr_glut
  key = bkey.decode("utf-8")
  if key == 'q':
    exit()
  glutPostRedisplay()

def special(key, x, y):
  global wmngr_glut
  wmngr_glut.special(key, x, y)

def mouse(button, state, x, y):
  global wmngr_glut
  wmngr_glut.mouse(button, state, x, y)

def motion(x, y):
  global wmngr_glut
  wmngr_glut.motion(x, y)

def main():
  global wmngr_glut

  ###################
  # GLUT Window Initialization
  glutInit()
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)  # zBuffer
  glutInitWindowSize(600, 600)
  glutInitWindowPosition(100, 100)
  glutCreateWindow("Simple GLUT")

  # Register callbacks
  glutReshapeFunc(reshape)
  glutDisplayFunc(display)
  glutMouseFunc(mouse)
  glutMotionFunc(motion)
  glutKeyboardFunc(keyboard)
  glutSpecialFunc(special)
  glutIdleFunc(idle)

  wmngr_glut = dfm2.glut.WindowManagerGLUT(0.3)

  dfm2.setSomeLighting()

  ####
  # Turn the flow of control over to GLUT
  glutMainLoop()

if __name__ == "__main__":
  main()
