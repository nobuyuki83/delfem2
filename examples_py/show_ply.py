from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

mshelm = dfm2.MeshElem("../test_inputs/bunny_2k.ply");
aabb = dfm2.AABB3( mshelm.minmax_xyz() )
dfm2.winDraw3d([mshelm,aabb])