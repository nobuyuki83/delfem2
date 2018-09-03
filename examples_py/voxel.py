from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

mshelm = None

def main():
  global mshelm

  voxelgrid = dfm2.VoxelGrid()
  voxelgrid.add(0, 0, 0)
  voxelgrid.add(1, 0, 0)
  voxelgrid.add(0, 1, 0)

  mshelm = dfm2.meshQuad3d_voxelGrid(voxelgrid)
  mshelm = mshelm.subdiv()
  mshelm = mshelm.subdiv()
  mshelm = mshelm.subdiv()

  axis = dfm2.AxisXYZ()

  dfm2.winDraw3d([mshelm,axis],(400,400))


if __name__ == "__main__":
  main()