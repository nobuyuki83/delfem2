from OpenGL.GL import *
from OpenGL.GLUT import *

import numpy as np

import sys
sys.path.append("../python")
import dfm2

mshelm = None

def display():
  glEnable(GL_LIGHTING)
  mshelm.draw()

  glDisable(GL_LIGHTING)
  glLineWidth(3)
  dfm2.gl.draw_axis(size=0.5)

def main():
  global mshelm
  mshelm = dfm2.dfm2.MeshElem("../test_inputs/bunny_2k.ply");
  mshelm.scale_xyz(0.01)

  win = dfm2.glut.GlutWindow(1.0)
  dfm2.dfm2.set_some_lighting()
  win.draw_loop(display)

if __name__ == "__main__":
  main()