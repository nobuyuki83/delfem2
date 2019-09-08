####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import numpy, random, pytest
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw
import PyDelFEM2.gl._gl

class Test_PBD():
  def test_pbd_tri(self,request):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D(edge_length=0.1)
    mesh = mesher.meshing(cad)
    pbd = dfm2.PBD()
    pbd.updated_topology(mesh)
    npIdP = cad.points_edge([0], mesh.np_pos)
    pbd.vec_bc[npIdP] = 1
    fvs = dfm2.FieldValueSetter("0.3*sin(2*t)", pbd.vec_val, 0,
                                mesh=mesh, npIdP=npIdP, dt=pbd.dt)
    mesh_def = dfm2.Mesh(np_pos=pbd.vec_val, np_elm=mesh.np_elm)
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([fvs,pbd,mesh_def],nframe=100)

  def test_pbd_hex(self,request):
    voxelgrid = dfm2.VoxelGrid()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(2, 0, 0)
    voxelgrid.add(1, 1, 0)
    msh = voxelgrid.mesh_hex()
    pbd = dfm2.PBD()
    pbd.updated_topology(msh)
    npIdP = numpy.array([0, 1, 2, 3], dtype=numpy.int32)
    pbd.vec_bc[npIdP] = 1
    fvs = dfm2.FieldValueSetter("0.4*sin(0.8*t)", pbd.vec_val, 1,
                                mesh=msh, npIdP=npIdP, dt=pbd.dt)
    ##
    mesh_def = dfm2.Mesh(np_pos=pbd.vec_val, np_elm=msh.np_elm, elem_type=dfm2.HEX)
    if request.config.getoption('--is_gl') == "true":                      
      dfm2.gl.glfw.winDraw3d([fvs,pbd,mesh_def],
                             nframe=100, camera_rotation=[0.1, 0.2, 0.0])

class Test_FEMPoission2D():
  def test1(self,request):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
    msh = mesher.meshing(cad)
    fem = dfm2.FEM_Poisson(source=1.0)
    fem.updated_topology(msh)
    npIdP = cad.points_edge([0,1,2,3], msh.np_pos)
    fem.ls.bc[npIdP] = 1
    fem.solve()
    ##
    field = dfm2.gl.VisFEM_ColorContour(fem, "value")
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([field],nframe=10)


class Test_FEMPoission3D():
  def test1(self,request):
    sdf0 = dfm2.CppSDF3_Sphere(0.55, [-0.5, 0, 0], True)
    sdf1 = dfm2.CppSDF3_Sphere(0.55, [+0.5, 0, 0], True)
    np_xyz, np_tet = dfm2.isosurface([sdf0,sdf1])
    msh = dfm2.Mesh(np_xyz, np_tet, dfm2.TET)
    npIdP0 = numpy.where(msh.np_pos[:,0]>+1)
    npIdP1 = numpy.where(msh.np_pos[:,0]<-1)
    fem = dfm2.FEM_Poisson()
    fem.updated_topology(msh)
    fem.ls.bc[npIdP0] = 1
    fem.ls.bc[npIdP1] = 2
    fem.value[:] = 0.5
    fem.value[npIdP0] = 0.0
    fem.value[npIdP1] = 1.0
    fem.solve()
    ##
    field = dfm2.gl.VisFEM_ColorContour(fem, "value")
    field.set_color_minmax()
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([field],nframe=10)


class Test_FEMDiffuse2D():
  def test1(self,request):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
    msh = mesher.meshing(cad)
    fem = dfm2.FEM_Diffuse()
    fem.updated_topology(msh)
    npIdP = cad.points_edge([0, 1, 2, 3], msh.np_pos)
    fem.ls.bc[npIdP] = 1
    ##
    field = dfm2.gl.VisFEM_ColorContour(fem, "value")
    field.color_min = 0.0
    field.color_max = 1.0
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([fem,field],nframe=100)


class Test_FemDiffuse3D():
  def test1(self,request):
    sdf0 = dfm2.CppSDF3_Sphere(0.55, [-0.5, 0, 0], True)
    sdf1 = dfm2.CppSDF3_Sphere(0.55, [+0.5, 0, 0], True)
    np_xyz, np_tet = dfm2.isosurface([sdf0,sdf1])
    msh = dfm2.Mesh(np_xyz, np_tet, dfm2.TET)
    npIdP0 = numpy.where(msh.np_pos[:,0]>+1)
    npIdP1 = numpy.where(msh.np_pos[:,0]<-1)
    fem = dfm2.FEM_Diffuse()
    fem.updated_topology(msh)
    fem.ls.bc[npIdP0] = 1
    fem.ls.bc[npIdP1] = 2
    fem.value[:] = 0.5
    fem.value[npIdP0] = 0.0
    fem.value[npIdP1] = 1.0
    ##
    field = dfm2.gl.VisFEM_ColorContour(fem, "value")
    field.color_min = 0.0
    field.color_max = 1.0
    if request.config.getoption('--is_gl') == "true":                      
      dfm2.gl.glfw.winDraw3d([fem,field],nframe=100)


class Test_FEMSolidLLinearStatic2D():
  def test1(self,request):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D(edge_length=0.1)
    msh = mesher.meshing(cad)
    fem = dfm2.FEM_SolidLinearStatic()
    fem.param_gravity_y = -0.1
    fem.updated_topology(msh)
    npIdP = cad.points_edge([3], msh.np_pos)
    fem.ls.bc[npIdP,:] = 1
    fem.solve()
    ##
    field = dfm2.gl.VisFEM_ColorContour(fem,name_disp="vec_val")
    if request.config.getoption('--is_gl') == "true":                      
      dfm2.gl.glfw.winDraw3d([field],nframe=10)


class Test_FEMSolidLLinearDynamic2D():
  def test1(self,request):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
    msh = mesher.meshing(cad)
    fem = dfm2.FEM_SolidLinearDynamic()
    fem.param_gravity_y = -0.1
    fem.updated_topology(msh)
    npIdP = cad.points_edge([3], msh.np_pos)
    fem.ls.bc[npIdP, :] = 1
    ##
    field = dfm2.gl.VisFEM_ColorContour(fem,name_disp="vec_val")
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([fem,field],nframe=100)


class Test_FEMSorkes2D():
  def test1(self,request):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
    msh = mesher.meshing(cad)
    fem = dfm2.FEM_StorksStatic2D(msh)
    npIdP0 = cad.points_edge([0,1,2,3], msh.np_pos)
    fem.ls.bc[npIdP0,0:2] = 1
    npIdP1 = cad.points_edge([2], msh.np_pos)
    fem.vec_val[npIdP1,0] = 1.0
    fem.solve()
    ####
    field_p = dfm2.gl.VisFEM_ColorContour(fem, name_color="vec_val", idim=2)
    field_p.set_color_minmax()
    field_v = dfm2.gl.VisFEM_Hedgehog(fem, name_vector="vec_val")
    axis = dfm2.gl.AxisXYZ(1.0)
    if request.config.getoption('--is_gl') == "true":                      
      dfm2.gl.glfw.winDraw3d([field_p, field_v, axis],nframe=10)


class Test_FEMCloth():
  def test1(self,request):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    
    mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
    msh = mesher.meshing(cad)
    fem = dfm2.FEM_Cloth()
    fem.updated_topology(msh)
    npIdP = cad.points_edge([2], msh.np_pos)
    fem.ls.bc[npIdP,0:3] = 1
    ####
    mesh2 = dfm2.Mesh(np_pos=fem.vec_val, np_elm=msh.np_elm)
    axis = dfm2.gl.AxisXYZ(1.0)
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([fem, mesh2, axis],nframe=100)

