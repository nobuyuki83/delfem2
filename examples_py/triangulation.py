from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

aXYVertexOuterLoop = [0,0, 1,0, 1,1, 0,1]
output = dfm2.triangulation(aXYVertexOuterLoop,edge_length=0.02)
mshelm = output.meshElem
dfm2.winDraw3d(mshelm)

