from OpenGL.GL import *

import sys
sys.path.append("../python")
import dfm2

if __name__ == "__main__":

  outerloop = [0,0, 1,0, 1,1, 0,1]
  output = dfm2.triangulation(outerloop,edge_length=0.02)
  mshelm = output.meshElem
  dfm2.winDraw3d(mshelm)

