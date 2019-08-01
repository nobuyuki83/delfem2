####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import PyDelFEM2 as dfm2
import PyDelFEM2.gl._glfw

import numpy


def example1():
  cad = dfm2.Cad2D()
  cad.add_polygon([+0,0, +1,0, +1,+1, 0,+1.0])
  cad.add_polygon([+2,0, +3,0, +3,+1, 2,+1.0])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  mesh = mesher.meshing(cad)

  ####
  pbd = dfm2.PBD_Cloth()
  pbd.param_gravity_y = -0.1
  pbd.dt = 0.08
  pbd.updated_topology(mesh)

  trans0 = dfm2.Trans_Rigid2DTo3D()
  trans0.org2 = numpy.array([0.5,0.5])
  trans0.org3 = numpy.array([0.0,0.0,0.5])

  trans1 = dfm2.Trans_Rigid2DTo3D()
  trans1.org2 = numpy.array([2.5,0.5])
  trans1.org3 = numpy.array([0.0,0.0,-0.5])
  trans1.R = dfm2.util.mat3_rot_cartesian(numpy.array([0,3.1415,0]))

  npIndP_Face0 = mesher.points_on_faces([0],cad)
  pbd.vec_val[npIndP_Face0] = trans0.trans(pbd.dmsh.np_pos[npIndP_Face0])

  npIndP_Face1 = mesher.points_on_faces([1],cad)
  pbd.vec_val[npIndP_Face1] = trans1.trans(pbd.dmsh.np_pos[npIndP_Face1])

  npIndP_Edge0a = mesher.points_on_one_edge(1,True,cad)
  npIndP_Edge0b = mesher.points_on_one_edge(7,True,cad)
  npIndP_Seam0 = numpy.vstack([npIndP_Edge0a,npIndP_Edge0b[::-1]]).transpose()
  pbd.elems_seam = npIndP_Seam0

  npIdP = mesher.points_on_edges([2,6],cad)
  pbd.bc[npIdP] = 1
  ####
  mesh2 = dfm2.Mesh(np_pos=pbd.vec_val,np_elm=mesh.np_elm)
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl._glfw.winDraw3d([pbd,mesh2,axis])


if __name__ == "__main__":
  example1()
