from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2
import dfm2.glfw

msh = dfm2.Mesh()
msh.read("../test_inputs/bunny_2k.ply")
aabb = dfm2.AABB3( msh.minmax_xyz() )
dfm2.glfw.winDraw3d([msh,aabb])

msh.read("../test_inputs/bunny_1k.obj");
aabb = dfm2.AABB3( msh.minmax_xyz() )
dfm2.glfw.winDraw3d([msh,aabb])