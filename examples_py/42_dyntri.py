####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import random
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def main_CppMeshDynTri2D_0():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1.0])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.1)
  dmesh = mesher.meshing(cad)
  dmesh.cdmsh.check()
  dfm2.gl.glfw.winDraw3d([dmesh.cdmsh,dfm2.gl.AxisXYZ()],winsize=(400,300))


def main_CppMeshDynTri2D_1():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1.0])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.1)
  dmesh = mesher.meshing(cad)
  cdmesh = dmesh.cdmsh
  cdmesh.check()
  for itr in range(100):
    itri0 = random.randint(0,cdmesh.ntri()-1)
    r0 = random.uniform(0.02, 0.98)
    r1 = random.uniform(0.01, 0.99-r0)
    ipo= cdmesh.insert_point_elem(itri0,r0,r1)
    cdmesh.delaunay_around_point(ipo)
    cdmesh.check()
  dfm2.gl.glfw.winDraw3d([cdmesh],winsize=(400,300))


def main_CppMeshDynTri3D():
  msh = dfm2.Mesh()
  msh.read("../test_inputs/bunny_2k.ply")
  dmesh = dfm2.CppMeshDynTri3D()
  dfm2.meshdyntri3d_initialize(dmesh, msh.np_pos, msh.np_elm)
  dmesh.check()
  for itr in range(500):
    itri0 = random.randint(0,dmesh.ntri()-1)
    iedge0 = random.randint(0,2)
    dmesh.delete_tri_edge(itri0,iedge0)
  dfm2.gl.glfw.winDraw3d([dmesh],winsize=(400,300))


def main_MeshDynTri2D_0():
  dmsh = dfm2.MeshDynTri2D()
  dmsh.meshing_loops([[-1,-1, +1,-1, +1,+1, -1,+1],
                      [-0.5, -0.5, +0.5, -0.5, 0.0, +0.5]],
                     edge_length=0.2)
  dfm2.gl.glfw.winDraw3d([dmsh])


def main_MeshDynTri2D_1():
  dmsh = dfm2.MeshDynTri2D()
  dmsh.meshing_loops([[-1,-1, +1,-1, +1,+1, -1,+1]],
                     edge_length=0.2)
  dfm2.gl.glfw.winDraw3d([dmsh])
  res = dmsh.refine_EdgeLongerThan_InsideCircle(0.05, 0.0,0.0, 0.5)
  dfm2.gl.glfw.winDraw3d([dmsh])


def main_MeshDynTri2D_2():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  dmsh = mesher.meshing(cad)
  ####
  fem = dfm2.FEM_ScalarPoisson(source=1.0)
  fem.updated_topology(dmsh)
  npIdP = cad.points_edge([0,1,2,3], dmsh.np_pos)
  fem.value[npIdP] = 0.0
  fem.ls.bc[npIdP] = 1
  fem.solve()
  vis = dfm2.gl.VisFEM_ColorContour(fem,name_color="value")
  dfm2.gl.glfw.winDraw3d([vis,dmsh])
  #####
  mapper = dmsh.refine_EdgeLongerThan_InsideCircle(0.05, 0.0,0.0,0.5)
  fem.updated_topology(dmsh,mapper=mapper)
  dfm2.gl.glfw.winDraw3d([vis,dmsh])
  #####
  fem.value[npIdP] = 0.0
  fem.ls.bc[npIdP] = 1
  fem.solve()
  dfm2.gl.glfw.winDraw3d([vis,dmsh])


def main_MeshDynTri2D_3():
  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  dmsh = mesher.meshing(cad)
  fem = dfm2.FEM_ShellCloth()
  fem.updated_topology(dmsh)
  npIdP = cad.points_edge([0], dmsh.np_pos)
  fem.ls.bc[npIdP,:] = 1
  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=dmsh.np_elm)
  dfm2.gl.glfw.winDraw3d([fem,mesh2])
  #####
  mesh2 = None
  mapper = dmsh.refine_EdgeLongerThan_InsideCircle(0.05, 0.0,0.0,0.5)
  fem.updated_topology(dmsh,mapper=mapper)
  fem.ls.bc[npIdP,:] = 1
  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=dmsh.np_elm)
  dfm2.gl.glfw.winDraw3d([fem,mesh2])
  '''  
  #####
  fem.vec_val[npIdP] = 0.0
  fem.ls.vec_bc[npIdP] = 1
  fem.solve()
  field = dfm2.Field(dmsh,val_color=fem.vec_val[:,0])
  dfm2.winDraw3d([field,dmsh])
  '''


if __name__ == "__main__":
  main_CppMeshDynTri2D_0()
  main_CppMeshDynTri2D_1()

  main_CppMeshDynTri3D()

  main_MeshDynTri2D_0()
  main_MeshDynTri2D_1()
  main_MeshDynTri2D_2()
  main_MeshDynTri2D_3()