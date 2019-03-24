from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

def draw_func():
  glEnable(GL_LIGHTING)
  msh.draw()

msh = dfm2.mesh_read("../test_inputs/bunny_2k.ply");
msh.scale_xyz(0.03)

win = dfm2.WindowGLFW(1.0,winsize=(400,300));
dfm2.setSomeLighting()
win.draw_loop([draw_func])
