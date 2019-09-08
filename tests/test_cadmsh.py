####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

#import unittest,
import pytest
import numpy, random, os, sys
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

os.chdir(os.path.dirname(os.path.abspath(__file__))) # for python3 setup.py test

class Test_CppCad2D(unittest.TestCase):
  def test1(self):
    ccad = dfm2.CppCad2D()
    ccad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
    self.assertTrue( ccad.check() )
    dfm2.gl.glfw.winDraw3d([ccad],nframe=20)

    ccad.add_vtx_face(0.0, 0.0, 0)
    self.assertTrue( ccad.check() )
    dfm2.gl.glfw.winDraw3d([ccad], nframe=20)

    ccad.add_vtx_edge(0.0,0.8, 2)
    self.assertTrue( ccad.check() )
    dfm2.gl.glfw.winDraw3d([ccad],nframe=20)

    ccad.set_edge_type(0, 1, [0.2, 0.3, -0.2, 0.3])
    self.assertTrue(ccad.check())
    dfm2.gl.glfw.winDraw3d([ccad], nframe=20)

    ccad.set_edge_type(0, 0, [])
    self.assertTrue(ccad.check())
    dfm2.gl.glfw.winDraw3d([ccad], nframe=20)


class Test_Cad2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
    self.assertEqual(cad.ccad.ind_vtx_face(0),[0,1,2,3])
    self.assertEqual(cad.ccad.ind_edge_face(0),[(0,True),(1,True),(2,True),(3,True)])
    self.assertTrue(cad.ccad.check())
    ####
    cad.add_vtx_edge(0,[0.0,0.0])
    self.assertTrue(cad.ccad.check())
    self.assertEqual(cad.ccad.ind_vtx_face(0),[0,4,1,2,3])
    self.assertEqual(cad.ccad.ind_edge_face(0),[(0,True),(4,True),(1,True),(2,True),(3,True)])

  def test2(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
    self.assertTrue(cad.ccad.check())
    ####
    cad.pick(-1,-1, 2.0)
    self.assertEqual( cad.ivtx_picked(), 0 )
    self.assertEqual( cad.iedge_picked(), -1 )
    self.assertEqual( cad.iface_picked(), -1 )
    ####
    cad.pick(+0,-1, 2.0)
    self.assertEqual( cad.ivtx_picked(), -1 )
    self.assertEqual( cad.iedge_picked(), 0 )
    self.assertEqual( cad.iface_picked(), -1 )
    ####
    cad.pick(+0.0,-0.1, 2.0)
    self.assertEqual( cad.ivtx_picked(), -1 )
    self.assertEqual( cad.iedge_picked(), -1 )
    self.assertEqual( cad.iface_picked(), 0 )



  def test3(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
    mesher = dfm2.Mesher_Cad2D()
    dmsh = mesher.meshing(cad)
#    msh,map_cad2msh = cad.mesh(0.02)
    self.assertEqual(dmsh.np_pos.shape[1],2)
    np_xy_bound = numpy.array(cad.ccad.xy_vtxctrl_face(0)).reshape([-1, 2])
    W = dfm2.mvc(dmsh.np_pos, np_xy_bound)
    self.assertEqual(W.ndim,2)
    self.assertEqual(W.shape[0],dmsh.np_pos.shape[0])
    self.assertEqual(W.shape[1],len(cad.ccad.xy_vtxctrl_face(0))/2)
    self.assertLess(numpy.linalg.norm(W.sum(axis=1)-numpy.ones((W.shape[0]))),1.0e-3)


class Test_MathExpression(unittest.TestCase):

  def test1(self):
    mee = dfm2.MathExpressionEvaluator()
    mee.set_key("x",3)
    mee.set_expression("x+3")
    mee.set_key("x",6)
    self.assertLess( mee.eval()-9.0, 1.0e-30 )


class Test_Mesh(unittest.TestCase):
  def __init__(self,str_method,is_gl=True):
    super().__init__(str_method)
    self.is_gl = is_gl

  def test1(self):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_2k.ply")
    self.assertIsNot(msh,None)
    if self.is_gl:
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test2(self):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_1k.obj")
    self.assertIsNot(msh,None)
    if self.is_gl:
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test3(self):
    msh = dfm2.Mesh()
    msh.set_grid((32,16))
    self.assertIsNot(msh,None)
    if self.is_gl:
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test4(self):
    voxelgrid = dfm2.VoxelGrid()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = voxelgrid.mesh_quad()
    msh = msh.subdiv()
    msh = msh.subdiv()
    msh = msh.subdiv()
    self.assertIsNot(msh,None)
    if self.is_gl:
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test5(self):
    voxelgrid = dfm2.VoxelGrid()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = voxelgrid.mesh_hex()
    msh = msh.subdiv()
    self.assertIsNot(msh,None)
    if self.is_gl:
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)


class Test_CppMeshDynTri3D(unittest.TestCase):
  def test0(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesher = dfm2.Mesher_Cad2D()
    mesh = mesher.meshing(cad)
    dmesh = dfm2.CppMeshDynTri3D()
    dfm2.meshdyntri3d_initialize(dmesh, mesh.np_pos, mesh.np_elm)
    dmesh.check()

  def test1(self):
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
    dfm2.gl.glfw.winDraw3d([dmesh], winsize=(400, 300), nframe=10)

  def test2(self):
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
    dfm2.gl.glfw.winDraw3d([dmesh], winsize=(400, 300), nframe=10)


class Test_MeshDynTri2D(unittest.TestCase):
  def test0(self):
    dmsh = dfm2.MeshDynTri2D()
    dmsh.meshing_loops([[-1,-1, +1,-1, +1,+1, -1,+1],
                        [-0.5, -0.5, +0.5, -0.5, 0.0, +0.5]],
                       edge_length=0.1)


class Test_GlTF(unittest.TestCase):
  def test0(self):
    gltf = PyDelFEM2.CppGLTF()
    # gltf.read("../test_inputs/RiggedFigure.glb")
    gltf.read("../test_inputs/CesiumMan.glb")
    # gltf.print()

    np_pos0, np_elm, np_rigw, np_rigj = dfm2.CppGLTF_GetMeshInfo(gltf, 0, 0)
    np_pos = np_pos0.copy()
    msh = dfm2.Mesh(np_pos, np_elm, dfm2.TRI)
    dfm2.gl.glfw.winDraw3d([msh], winsize=(400, 300), nframe=10)

    bone_array = dfm2.CppGLTF_GetBones(gltf, 0)
    bone_array.set_rotation_bryant(0, [-3.1415 * 0.5, 0.0, 0.0])
    bone_array.set_translation(0, [0.0, 0.0, +0.2])
    dfm2.update_rig_skin(np_pos,
                         np_pos0, np_elm, bone_array, np_rigw, np_rigj)
    msh = dfm2.Mesh(np_pos, np_elm, dfm2.TRI)
    axis = dfm2.gl.AxisXYZ(1)
    dfm2.gl.glfw.winDraw3d([msh, axis], winsize=(400, 300), nframe=10)