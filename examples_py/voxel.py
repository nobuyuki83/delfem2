from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

def show_hex(voxelgrid):
  msh = voxelgrid.mesh_hex3d()
  msh = msh.subdiv()
  axis = dfm2.AxisXYZ()
  dfm2.winDraw3d([msh, axis], (400, 400))

def show_quad(voxelgrid):
  msh = voxelgrid.mesh_quad3d()
  msh = msh.subdiv()
  msh = msh.subdiv()
  msh = msh.subdiv()
  axis = dfm2.AxisXYZ()
  dfm2.winDraw3d([msh, axis], (400, 400))

def main():
  voxelgrid = dfm2.Grid3D()
  voxelgrid.add(0, 0, 0)
  voxelgrid.add(1, 0, 0)
  voxelgrid.add(0, 1, 0)
  ###
  show_quad(voxelgrid)
  show_hex(voxelgrid)


if __name__ == "__main__":
  main()