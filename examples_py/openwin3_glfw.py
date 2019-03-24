from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

msh = dfm2.mesh_read("../test_inputs/bunny_2k.ply");
msh.scale_xyz(0.03)
dfm2.winDraw3d([msh],winsize=(400,300))