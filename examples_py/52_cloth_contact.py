####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import sys
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw

def example1():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh,map_cad2mesh = cad.mesh(0.05)
  fem = dfm2.FEM_Cloth()
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([2], mesh.np_pos)
  fem.ls.bc[npIdP,0:3] = 1
  fem.sdf.list_sdf.append( dfm2.SDF_Sphere(0.55,[0,+0.5,-1.0],True) )
  ####
  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=mesh.np_elm)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.glfw.winDraw3d([fem,mesh2,axis,fem.sdf])


if __name__ == "__main__":
  example1()
