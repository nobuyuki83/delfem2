####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import OpenGL.GL as gl
import OpenGL.GLUT as glut

import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glut

wmngr_glut = None

def display():
  global wmngr_glut
  gl.glClearColor(0.3, 0.5, 0.8, 1.0)
  gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
  gl.glEnable(gl.GL_DEPTH_TEST)
  
  wmngr_glut.camera.set_gl_camera()

  gl.glEnable(gl.GL_LIGHTING)
  glut.glutSolidTeapot(0.1)

  gl.glDisable(gl.GL_LIGHTING)
  gl.glLineWidth(3)
  glut.glutSwapBuffers()

def reshape(width, height):
  gl.glViewport(0, 0, width, height)

def idle():
  glut.glutPostRedisplay()

def keyboard(bkey, x, y):
  global wmngr_glut
  key = bkey.decode("utf-8")
  if key == 'q':
    exit()
  glut.glutPostRedisplay()

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
  glut.glutInit()
  glut.glutInitDisplayMode(glut.GLUT_DOUBLE | glut.GLUT_RGB | glut.GLUT_DEPTH)  # zBuffer
  glut.glutInitWindowSize(600, 600)
  glut.glutInitWindowPosition(100, 100)
  glut.glutCreateWindow("Simple GLUT")

  # Register callbacks
  glut.glutReshapeFunc(reshape)
  glut.glutDisplayFunc(display)
  glut.glutMouseFunc(mouse)
  glut.glutMotionFunc(motion)
  glut.glutKeyboardFunc(keyboard)
  glut.glutSpecialFunc(special)
  glut.glutIdleFunc(idle)

  wmngr_glut = dfm2.gl.glut.WindowManagerGLUT(0.3)

  dfm2.gl.setSomeLighting()

  ####
  # Turn the flow of control over to GLUT
  glut.glutMainLoop()

if __name__ == "__main__":
  main()
