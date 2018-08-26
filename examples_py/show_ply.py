from OpenGL.GL import *

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
  mshelm = dfm2.dfm2.MeshElem("../test_inputs/bunny_2k.ply");
  mshelm.scaleXYZ(0.01)

  win = dfm2.gl_fw.WindowGLFW(0.5)
  dfm2.dfm2.setSomeLighting()
  win.draw_loop(display)

if __name__ == "__main__":
  main()