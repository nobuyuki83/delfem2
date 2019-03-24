from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

msh = None

def main():
  global msh

  voxelgrid = dfm2.VoxelGrid()
  voxelgrid.add(0, 0, 0)
  voxelgrid.add(1, 0, 0)
  voxelgrid.add(0, 1, 0)

  msh = dfm2.mesh_voxelgrid(voxelgrid=voxelgrid)
  msh = msh.subdiv()
  msh = msh.subdiv()
  msh = msh.subdiv()

  axis = dfm2.AxisXYZ()

  dfm2.winDraw3d([msh, axis], (400, 400))


if __name__ == "__main__":
  main()