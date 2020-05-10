####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import pytest
import numpy, random, os, sys, math
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

os.chdir(os.path.dirname(os.path.abspath(__file__))) # for python3 setup.py test


class Test_MathExpression():
  def test1(self):
    mee = dfm2.MathExpressionEvaluator()
    mee.set_key("x",3)
    mee.set_expression("x+3")
    mee.set_key("x",6)
    assert math.fabs(mee.eval()-9.0)<1.0e-30


class Test_Mesh():

  def test1(self,request):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_2k.ply")
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test2(self,request):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_1k.obj")
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test3(self,request):
    msh = dfm2.Mesh()
    msh.set_grid((32,16))
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test4(self,request):
    voxelgrid = dfm2.VoxelGrid()
    voxelgrid.initialize(2,2,1)
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = voxelgrid.mesh_quad()
    msh = msh.subdiv()
    msh = msh.subdiv()
    msh = msh.subdiv()
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":    
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test5(self,request):
    voxelgrid = dfm2.VoxelGrid()
    voxelgrid.initialize(2,2,1)
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = voxelgrid.mesh_hex()
    msh = msh.subdiv()
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":            
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)



class Test_CppMeshDynTri3D():
  def test0(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D()
    mesh = mesher.meshing(cad)
    cdmsh = dfm2.CppMeshDynTri3D()
    dfm2.meshdyntri3d_initialize(cdmsh, mesh.np_pos, mesh.np_elm)
    cdmsh.check()

  def test1(self,request):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D(edge_length=0.08)
    mesh = mesher.meshing(cad)
    dmesh = dfm2.CppMeshDynTri3D()
    dfm2.meshdyntri3d_initialize(dmesh, mesh.np_pos, mesh.np_elm)
    dmesh.check()
    for itr in range(100):
      itri0 = random.randint(0, dmesh.ntri() - 1)
      r0 = random.uniform(0.02, 0.98)
      r1 = random.uniform(0.01, 0.99 - r0)
      ipo = dmesh.insert_point_elem(itri0, r0, r1)
      dmesh.delaunay_around_point(ipo)
      dmesh.check()
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([dmesh], winsize=(400, 300), nframe=10)

  def test2(self,request):
    print(os.path.abspath(__file__))
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_2k.ply")
    dmesh = dfm2.CppMeshDynTri3D()
    dfm2.meshdyntri3d_initialize(dmesh, msh.np_pos, msh.np_elm)
    dmesh.check()
    for itr in range(100):
      itri0 = random.randint(0, dmesh.ntri() - 1)
      iedge0 = random.randint(0, 2)
      dmesh.delete_tri_edge(itri0, iedge0)
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([dmesh], winsize=(400, 300), nframe=10)


class Test_MeshDynTri2D():
  def test0(self):
    dmsh = dfm2.MeshDynTri2D()
    dmsh.meshing_loops([[-1,-1, +1,-1, +1,+1, -1,+1],
                        [-0.5, -0.5, +0.5, -0.5, 0.0, +0.5]],
                       edge_length=0.1)




class Test_GlTF():
  def test0(self,request):
    gltf = PyDelFEM2.CppGLTF()
    # gltf.read("../test_inputs/RiggedFigure.glb")
    gltf.read("../test_inputs/CesiumMan.glb")
    # gltf.print()

    np_pos0, np_elm, np_rigw, np_rigj = dfm2.CppGLTF_GetMeshInfo(gltf, 0, 0)
    np_pos = np_pos0.copy()
    msh = dfm2.Mesh(np_pos, np_elm, dfm2.TRI)
    if request.config.getoption('--is_gl') == "true":                      
      dfm2.gl.glfw.winDraw3d([msh], winsize=(400, 300), nframe=10)

    bone_array = dfm2.CppGLTF_GetBones(gltf, 0)
    bone_array.set_rotation_bryant(0, [-3.1415 * 0.5, 0.0, 0.0])
    bone_array.set_translation(0, [0.0, 0.0, +0.2])
    dfm2.update_rig_skin(np_pos,
                         np_pos0, np_elm, bone_array, np_rigw, np_rigj)
    msh = dfm2.Mesh(np_pos, np_elm, dfm2.TRI)
    axis = dfm2.gl.AxisXYZ(1)
    if request.config.getoption('--is_gl') == "true":                  
      dfm2.gl.glfw.winDraw3d([msh, axis], winsize=(400, 300), nframe=10)


class Test_Tex():
  def test0(self,request):
    path_img = "../test_inputs/lenna.png"
    np_img = dfm2.imread(path_img)
    tex = dfm2.gl.get_texture(np_img,"rgb")
    print(tex.minmax_xyz())
    axis = dfm2.gl.AxisXYZ(100)
    if request.config.getoption('--is_gl') == "true":
      dfm2.gl.glfw.winDraw3d([tex, axis], winsize=(400, 300), nframe=10)
