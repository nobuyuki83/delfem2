####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import numpy
import sys
sys.path.append("..")
import pydelfem2 as dfm2
import pydelfem2.gl._glfw

def poisson(cad,mesh,map_cad2mesh):
  fem = dfm2.FEM_Poisson(source=1.0)
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.bc[npIdP] = 1
  fem.solve()
  field = dfm2.VisFEM_ColorContour(fem,"value")
  dfm2.gl._glfw.winDraw3d([field])
  field.write_vtk("poisson2d.vtk")


def poisson_ms(cad, mesh):
  npIdP_ms = cad.points_edge([0], mesh.np_pos)
  vec_ms = numpy.zeros((mesh.np_pos.shape[0],1),dtype=numpy.int32)
  vec_ms[:] = -1
  vec_ms[npIdP_ms] = npIdP_ms[0]
  vec_ms[npIdP_ms[0]] = -1
  fem = dfm2.FEM_Poisson(source=0.2)
  fem.updated_topology(mesh=mesh, master_slave_pattern=vec_ms)
  npIdP_fix = cad.points_edge([2], mesh.np_pos)
  fem.ls.bc[npIdP_fix] = 1
  fem.solve()
  field = dfm2.VisFEM_ColorContour(fem,"value")
  dfm2.gl._glfw.winDraw3d([field])


def diffuse(cad,mesh):
  fem = dfm2.FEM_Diffuse()
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([0,1,2,3], mesh.np_pos)

  fem.ls.bc[npIdP] = 1
  ####
  field = dfm2.VisFEM_ColorContour(fem,name_color="value")
  field.draw_val_min = 0.0
  field.draw_val_max = 0.3
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl._glfw.winDraw3d([fem,field,axis])


def linear_solid_static(cad,mesh):
  fem = dfm2.FEM_SolidLinearStatic()
  fem.param_gravity_y = -0.1
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([3], mesh.np_pos)
  fem.ls.bc[npIdP,:] = 1
  fem.solve()
  field = dfm2.VisFEM_ColorContour(fem,name_disp="vec_val")
  dfm2.gl._glfw.winDraw3d([field])
  field.write_vtk("linearsolid2d.vtk")


def linear_solid_eigen(mesh):
  fem = dfm2.FEM_SolidLinearEigen()
  fem.updated_topology(mesh)
  fem.ls.f[:] = numpy.random.uniform(-1,1, mesh.np_pos.shape )
  field = dfm2.VisFEM_ColorContour(fem,name_disp="mode")
  dfm2.gl._glfw.winDraw3d([fem,field])


def linear_solid_dynamic(cad,mesh):
  fem = dfm2.FEM_SolidLinearDynamic()
  fem.param_gravity_y = -0.1
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([3], mesh.np_pos)
  fem.ls.bc[npIdP,:] = 1
  print(fem.ls.conv_hist)
  ####
  field = dfm2.VisFEM_ColorContour(fem,name_disp="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl._glfw.winDraw3d([fem,field,axis])


def storks_static(cad,mesh):
  fem = dfm2.FEM_StorksStatic2D(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field_p = dfm2.VisFEM_ColorContour(fem, name_color="vec_val",idim=2)
  field_p.set_color_minmax()
  field_v = dfm2.VisFEM_Hedgehog(fem, name_vector="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl._glfw.winDraw3d([field_p,field_v,axis])

def storks_dynamic(cad,mesh):
  fem = dfm2.FEM_StorksDynamic2D(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field_p = dfm2.VisFEM_ColorContour(fem, name_color="vec_val",idim=2)
  field_p.set_color_minmax()
  field_v = dfm2.VisFEM_Hedgehog(fem, name_vector="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl._glfw.winDraw3d([fem,field_p,field_v,axis])


def navir_storks(cad,mesh):
  fem = dfm2.FEM_NavierStorks2D(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  ####
  field_p = dfm2.VisFEM_ColorContour(fem, name_color="vec_val",idim=2)
  field_p.set_color_minmax()
  field_v = dfm2.VisFEM_Hedgehog(fem, name_vector="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl._glfw.winDraw3d([fem,field_p,field_v,axis])


def fem_cloth():
  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, +0.8,+1, -0.8,+1, -1,+1])
  mesh,map_cad2mesh = cad.mesh(0.1)
  ####
  fem = dfm2.FEM_Cloth()
  fem.dt = 0.08
  fem.lmd = 1000
  fem.myu = 100
  fem.gravity = (0,-1,0.01)
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([2,4], mesh.np_pos)
  fem.ls.bc[npIdP,0:3] = 1
  ####
  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=mesh.np_elm)
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl._glfw.winDraw3d([fem,mesh2,axis])


def pbd1(cad,mesh):
  pbd = dfm2.PBD()
  pbd.updated_topology(mesh)
  npIdP = cad.points_edge([0], mesh.np_pos)
  pbd.vec_bc[npIdP] = 1
  fvs = dfm2.FieldValueSetter("0.3*sin(2*t)", pbd.vec_val, 0,
                              mesh=mesh, npIdP=npIdP, dt=pbd.dt)
  ####
  mesh2 = dfm2.Mesh(np_pos=pbd.vec_val,np_elm=mesh.np_elm)
  dfm2.gl._glfw.winDraw3d([fvs,pbd,mesh2])


def pbd_cloth():
  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, +0.8,+1, -0.8,+1, -1,+1])
  mesh,map_cad2mesh = cad.mesh(0.1)
  ####
  pbd = dfm2.PBD_Cloth()
  pbd.param_gravity_y = -0.1
  pbd.param_gravity_z = -0.001
  pbd.dt = 0.08
  pbd.updated_topology(mesh)
  npIdP = cad.points_edge([2,4], mesh.np_pos)
  pbd.bc[npIdP] = 1
  ####
  mesh2 = dfm2.Mesh(np_pos=pbd.vec_val,np_elm=mesh.np_elm)
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl._glfw.winDraw3d([pbd,mesh2,axis])


def main():
  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  mesh,map_cad2mesh = cad.mesh(0.05)
  poisson(cad,mesh,map_cad2mesh)
  diffuse(cad,mesh)
  linear_solid_static(cad,mesh)
  linear_solid_dynamic(cad,mesh)
  storks_static(cad,mesh)
  storks_dynamic(cad,mesh)
  navir_storks(cad,mesh)

  fem_cloth()

  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,0, +0,+0, 0,+1, -1,+1.0])
  mesh,map_cad2mesh = cad.mesh(0.05)
  poisson_ms(cad, mesh)

  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh,map_cad2mesh = cad.mesh(0.2)
  pbd1(cad,mesh)

  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-0.2, +1,-0.2, +1,+0.2, -1,+0.2])
  msh2,map_cad2mesh = cad.mesh(0.05)
  msh25 = dfm2.Mesh()
  msh25.set_extrude(msh2,1)
  msh25.np_pos[:,2] *= 0.05
  linear_solid_eigen(msh25)

  pbd_cloth()

if __name__ == "__main__":
  main()
