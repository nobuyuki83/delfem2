from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

cad = dfm2.Cad2D()
cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
dfm2.winDraw3d([cad])