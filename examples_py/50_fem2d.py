####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import numpy
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def femScalarPoisson(cad, mesh, mesher):
  fem = dfm2.FEM_ScalarPoisson(source=1.0)
  fem.updated_topology(mesh)
  npIdP = mesher.points_on_edges([0,1,2,3], cad)
  fem.ls.bc[npIdP] = 1
  fem.solve()
  field = dfm2.gl.VisFEM_ColorContour(fem,"value")
  dfm2.gl.glfw.winDraw3d([field])
  field.write_vtk("poisson2d.vtk")


def femScalarPoissonMasterSlave(cad, mesh):
  npIdP_ms = cad.points_edge([0], mesh.np_pos)
  vec_ms = numpy.zeros((mesh.np_pos.shape[0],1),dtype=numpy.int32)
  vec_ms[:] = -1
  vec_ms[npIdP_ms] = npIdP_ms[0]
  vec_ms[npIdP_ms[0]] = -1
  fem = dfm2.FEM_ScalarPoisson(source=0.2)
  fem.updated_topology(mesh=mesh, master_slave_pattern=vec_ms)
  npIdP_fix = cad.points_edge([2], mesh.np_pos)
  fem.ls.bc[npIdP_fix] = 1
  fem.solve()
  field = dfm2.gl.VisFEM_ColorContour(fem,"value")
  dfm2.gl.glfw.winDraw3d([field])


def femSclarDiffuse(cad, mesh, mesher):
  fem = dfm2.FEM_ScalarDiffuse()
  fem.updated_topology(mesh)
  npIdP = mesher.points_on_edges([0,1,2,3], cad)

  fem.ls.bc[npIdP] = 1
  ####
  field = dfm2.gl.VisFEM_ColorContour(fem,name_color="value")
  field.draw_val_min = 0.0
  field.draw_val_max = 0.3
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([fem,field,axis])

####################################################

def femSolidLinear_Static(cad, mesh):
  fem = dfm2.FEM_SolidLinearStatic()
  fem.param_gravity_y = -0.1
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([3], mesh.np_pos)
  fem.ls.bc[npIdP,:] = 1
  fem.solve()
  field = dfm2.gl.VisFEM_ColorContour(fem,name_disp="vec_val")
  dfm2.gl.glfw.winDraw3d([field])
  field.write_vtk("linearsolid2d.vtk")


def femSolidLinear_Eigen(mesh):
  fem = dfm2.FEM_SolidLinearEigen()
  fem.updated_topology(mesh)
  fem.ls.f[:] = numpy.random.uniform(-1,1, mesh.np_pos.shape )
  field = dfm2.gl.VisFEM_ColorContour(fem,name_disp="mode")
  dfm2.gl.glfw.winDraw3d([fem,field])


def femSolidLinear_Dynamic(cad, mesh):
  fem = dfm2.FEM_SolidLinearDynamic()
  fem.param_gravity_y = -0.1
  fem.updated_topology(mesh)
  npIdP = cad.points_edge([3], mesh.np_pos)
  fem.ls.bc[npIdP,:] = 1
  print(fem.ls.conv_hist)
  ####
  field = dfm2.gl.VisFEM_ColorContour(fem,name_disp="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([fem,field,axis])


############################

def femShellPlateBendingMitc3_Static(cad:dfm2.Cad2D, mesher, mesh:dfm2.Mesh):
  fem = dfm2.FEM_ShellPlateBendingMITC3()
  fem.param_gravity_z = -0.01
  fem.updated_topology(mesh)
  fem.ls.bc[mesher.points_on_one_edge(3,True,cad),:] = 1
  fem.solve()
  print(fem.ls.conv_hist,len(fem.ls.conv_hist))
  pos1 = numpy.zeros((mesh.np_pos.shape[0],3))
  pos1[:,:2] = mesh.np_pos
  pos1[:,2] = fem.disp[:,0]
  mesh1 = dfm2.Mesh(pos1,mesh.np_elm, dfm2.TRI)
  dfm2.gl.glfw.winDraw3d([mesh1],camera_rotation=[-1.2,0,0])


def femShellPlateBendingMitc3_Eigen(mesh):
  fem = dfm2.FEM_ShellPlateBendingMITC3_Eigen()
  fem.updated_topology(mesh)
  fem.ls.f[:] = numpy.random.uniform(-1,1, (mesh.np_pos.shape[0],3) )
  for itr in range(60):
    fem.solve()
  pos1 = numpy.zeros((mesh.np_pos.shape[0],3))
  pos1[:,:2] = mesh.np_pos
  pos1[:,2] = fem.mode[:,0]
  mesh1 = dfm2.Mesh(pos1,mesh.np_elm, dfm2.TRI)
  dfm2.gl.glfw.winDraw3d([mesh,mesh1],camera_rotation=[-1.2,0,0])

def femShellCloth():
  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, +0.8,+1, -0.8,+1, -1,+1])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  mesh = mesher.meshing(cad)
  ####
  fem = dfm2.FEM_ShellCloth()
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
  dfm2.gl.glfw.winDraw3d([fem,mesh2,axis])

##########################################

def femFluidStokes_Static(cad, mesh):
  fem = dfm2.FEM_FluidStorksStatic(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field_p = dfm2.gl.VisFEM_ColorContour(fem, name_color="vec_val",idim=2)
  field_p.set_color_minmax()
  field_v = dfm2.gl.VisFEM_Hedgehog(fem, name_vector="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([field_p,field_v,axis])


def femFluidStokes_Dynamic(cad, mesh):
  fem = dfm2.FEM_FluidStorksDynamic(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field_p = dfm2.gl.VisFEM_ColorContour(fem, name_color="vec_val",idim=2)
  field_p.set_color_minmax()
  field_v = dfm2.gl.VisFEM_Hedgehog(fem, name_vector="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([fem,field_p,field_v,axis])


def femFluidNavirStokes(cad, mesh):
  fem = dfm2.FEM_FluidNavierStorks(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  ####
  field_p = dfm2.gl.VisFEM_ColorContour(fem, name_color="vec_val",idim=2)
  field_p.set_color_minmax()
  field_v = dfm2.gl.VisFEM_Hedgehog(fem, name_vector="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([fem,field_p,field_v,axis])

#######################################################

def pbd1(cad,mesh):
  pbd = dfm2.PBD()
  pbd.updated_topology(mesh)
  npIdP = cad.points_edge([0], mesh.np_pos)
  pbd.vec_bc[npIdP] = 1
  fvs = dfm2.FieldValueSetter("0.3*sin(2*t)", pbd.vec_val, 0,
                              mesh=mesh, npIdP=npIdP, dt=pbd.dt)
  ####
  mesh2 = dfm2.Mesh(np_pos=pbd.vec_val,np_elm=mesh.np_elm)
  dfm2.gl.glfw.winDraw3d([fvs,pbd,mesh2])


def pbd_cloth():
  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, +0.8,+1, -0.8,+1, -1,+1])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  mesh = mesher.meshing(cad)
  ####
  pbd = dfm2.PBD_Cloth()
  pbd.param_gravity_y = -0.1
  pbd.param_gravity_z = -0.001
  pbd.dt = 0.08
  pbd.updated_topology(mesh)
  npIdP = mesher.points_on_edges([2,4], cad)
  pbd.bc[npIdP] = 1
  ####
  mesh2 = dfm2.Mesh(np_pos=pbd.vec_val,np_elm=mesh.np_elm)
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([pbd,mesh2,axis])


def main():

  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  mesh = mesher.meshing(cad)
  femScalarPoisson(cad, mesh, mesher)
  femSclarDiffuse(cad, mesh, mesher)
  femSolidLinear_Static(cad, mesh)
  femSolidLinear_Dynamic(cad, mesh)
  femFluidStokes_Static(cad, mesh)
  femFluidStokes_Dynamic(cad, mesh)
  femFluidNavirStokes(cad, mesh)

  femShellCloth()

  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,0, +0,+0, 0,+1, -1,+1.0])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  mesh = mesher.meshing(cad)
  femScalarPoissonMasterSlave(cad, mesh)

  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.1)
  mesh = mesher.meshing(cad)
  pbd1(cad,mesh)

  pbd_cloth()

  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1, -0.2, +1, -0.2, +1, +0.2, -1, +0.2])
  mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
  msh2 = mesher.meshing(cad)
  femShellPlateBendingMitc3_Static(cad, mesher, msh2)
  femShellPlateBendingMitc3_Eigen(msh2)
  #  return

  msh25 = dfm2.Mesh()
  msh25.set_extrude(msh2, 1)
  msh25.np_pos[:, 2] *= 0.05
  femSolidLinear_Eigen(msh25)

if __name__ == "__main__":
  main()
