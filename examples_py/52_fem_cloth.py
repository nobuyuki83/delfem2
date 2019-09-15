####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def example1():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  mesh = mesher.meshing(cad)

  fem = dfm2.FEM_ShellCloth()
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([2], mesh.np_pos)
  fem.ls.bc[npIdP,0:3] = 1
  fem.sdf = dfm2.CppSDF3_Sphere(0.55,[0,+0.5,-1.0],True)

  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=mesh.np_elm)
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([fem,mesh2,axis,fem.sdf])


if __name__ == "__main__":
  example1()
