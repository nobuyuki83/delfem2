from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

msh = dfm2.mesh_read("../test_inputs/bunny_2k.ply");
aabb = dfm2.AABB3( msh.minmax_xyz() )
dfm2.winDraw3d([msh,aabb])

msh = dfm2.mesh_read("../test_inputs/bunny_1k.obj");
aabb = dfm2.AABB3( msh.minmax_xyz() )
dfm2.winDraw3d([msh,aabb])