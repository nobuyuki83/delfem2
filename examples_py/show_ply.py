from OpenGL.GL import *

import sys
sys.path.append("../python")
import dfm2

mshelm = dfm2.MeshElem("../test_inputs/bunny_2k.ply");
mshelm.scaleXYZ(0.01)
dfm2.winDraw3d(mshelm)