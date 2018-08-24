from OpenGL.GL import *
from OpenGL.GLUT import *

import numpy as np

import sys
sys.path.append("../python")
import dfm2

mshelm = None

def display():
  glEnable(GL_LIGHTING)
  mshelm.drawFace_elemWiseNorm()

  glDisable(GL_LIGHTING)
  glLineWidth(3)  
  dfm2.gl.draw_axis(size=0.5)

def main():
  global mshelm

  voxelgrid = dfm2.dfm2.VoxelGrid()
  voxelgrid.add(0,0,0)
  voxelgrid.add(1,0,0)

  mshelm = dfm2.dfm2.meshQuad3d_voxelGrid(voxelgrid)

  win = dfm2.glut.GlutWindow(2.0)
  dfm2.dfm2.setSomeLighting()
  win.draw_loop(display)

if __name__ == "__main__":
  main()