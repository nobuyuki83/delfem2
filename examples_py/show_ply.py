from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

msh = dfm2.Mesh("../test_inputs/bunny_2k.ply");
aabb = dfm2.AABB3( msh.minmax_xyz() )
dfm2.winDraw3d([msh,aabb])