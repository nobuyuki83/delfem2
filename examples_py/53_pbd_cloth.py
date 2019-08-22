####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

import numpy

########################################

def example1():
  cad = dfm2.Cad2D()
  cad.add_polygon([+0,0, +1,0, +1,+1, 0,+1.0])
  cad.add_polygon([+2,0, +3,0, +3,+1, 2,+1.0])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  mesh = mesher.meshing(cad)

  ####
  pbd = dfm2.PBD_Cloth()
#  pbd.param_gravity_y = -0.1
  pbd.dt = 0.08
  pbd.updated_topology(mesh)
  pbd.sdf = dfm2.SDF()
  pbd.sdf.add( dfm2.CppSDF3_Sphere(0.3, [0.0, 0.0, 0.0], True) )

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

  npIndP_Edge1a = mesher.points_on_one_edge(3,True,cad)
  npIndP_Edge1b = mesher.points_on_one_edge(5,True,cad)
  npIndP_Seam1 = numpy.vstack([npIndP_Edge1a,npIndP_Edge1b[::-1]]).transpose()

  pbd.elems_seam = numpy.vstack([npIndP_Seam0,npIndP_Seam1]).copy().astype(numpy.uint32) # to allign data

  mesh2 = dfm2.Mesh(np_pos=pbd.vec_val,np_elm=mesh.np_elm)
  mesh3 = dfm2.Mesh(np_pos=pbd.vec_val, np_elm=pbd.elems_seam, elem_type=dfm2.LINE)

  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([pbd,pbd.sdf,mesh2,mesh3,axis])


  msh_trg = dfm2.Mesh()
  msh_trg.set_sphere(0.3, 16, 16)

  pbd.sdf = dfm2.Collider_PointsToMeshTri3D()
  pbd.sdf.set_mesh(msh_trg)

  dfm2.gl.glfw.winDraw3d([pbd,pbd.sdf,mesh2,mesh3,msh_trg,axis])



if __name__ == "__main__":
  example1()
