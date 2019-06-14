####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import sys, numpy
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw

def show_hex(voxelgrid):
  msh = voxelgrid.mesh_hex3d()
  msh = msh.subdiv()
  axis = dfm2.AxisXYZ()
  dfm2.glfw.winDraw3d([msh, axis], (400, 400))

def pbd_hex(voxelgrid):
  msh = voxelgrid.mesh_hex3d()
  pbd = dfm2.PBD(msh)
  ####
  npIdP = numpy.array([0,1,2,3],dtype=numpy.int32)
  pbd.vec_bc[npIdP] = 1
  fvs = dfm2.FieldValueSetter("0.4*sin(0.8*t)", pbd.vec_val, 1,
                              mesh=msh, npIdP=npIdP, dt=pbd.dt)
  ####
  msh_def = dfm2.Mesh(np_pos=pbd.vec_val,np_elm=msh.np_elm)
  axis = dfm2.AxisXYZ()
  dfm2.glfw.winDraw3d([fvs, pbd, msh_def, axis], (400, 400))




def show_quad(voxelgrid):
  msh = voxelgrid.mesh_quad3d()
  msh = msh.subdiv()
  msh = msh.subdiv()
  msh = msh.subdiv()
  axis = dfm2.AxisXYZ()
  dfm2.glfw.winDraw3d([msh, axis], (400, 400))

def main():
  grid3d = dfm2.Grid3D()
  grid3d.add(0, 0, 0)
  grid3d.add(1, 0, 0)
  grid3d.add(2, 0, 0)
  grid3d.add(1, 1, 0)
  ###
  show_quad(grid3d)
  show_hex(grid3d)
  pbd_hex(grid3d)

if __name__ == "__main__":
  main()