from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

mshelm = dfm2.MeshElem("../test_inputs/bunny_2k.ply");
mshelm.scaleXYZ(0.03)
dfm2.winDraw3d(mshelm,winsize=(400,300))