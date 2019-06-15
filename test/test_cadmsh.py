import unittest, numpy, random
import sys
sys.path.append("../module_py")
import delfem2 as dfm2

class Test_Cad2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
    cad.ccad.check()
    cad.add_point_edge(0.0,0.0, 0)

  def test2(self):
    cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])

  def test3(self):
    cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
    msh = cad.mesh(0.02)
    self.assertEqual(msh.np_pos.shape[1],2)
    W = dfm2.WeightMVC_CadMesh(cad.ccad, msh)
    self.assertEqual(W.ndim,2)
    self.assertEqual(W.shape[0],msh.np_pos.shape[0])
    self.assertEqual(W.shape[1],len(cad.ccad.getVertexXY_face(0))/2)
    self.assertLess(numpy.linalg.norm(W.sum(axis=1)-numpy.ones((W.shape[0]))),1.0e-3)


class Test_MathExpression(unittest.TestCase):
  def test1(self):
    mee = dfm2.MathExpressionEvaluator()
    mee.set_key("x",3)
    mee.set_expression("x+3")
    mee.set_key("x",6)
    self.assertLess( mee.eval()-9.0, 1.0e-30 )


class Test_Mesh(unittest.TestCase):
  def test1(self):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_2k.ply")
    self.assertIsNot(msh,None)

  def test2(self):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_1k.obj");
    self.assertIsNot(msh,None)

  def test3(self):
    msh = dfm2.mesh_grid((32,64))
    self.assertIsNot(msh,None)

  def test4(self):
    voxelgrid = dfm2.Grid3D()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = voxelgrid.mesh_quad3d()
    msh = msh.subdiv()
    msh = msh.subdiv()
    msh = msh.subdiv()
    self.assertIsNot(msh,None)

  def test5(self):
    voxelgrid = dfm2.Grid3D()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = voxelgrid.mesh_hex3d()
    msh = msh.subdiv()
    self.assertIsNot(msh,None)

class Test_CppMeshDynTri3D(unittest.TestCase):
  def test0(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1.0])
    mesh = cad.mesh(0.1)
    dmesh = dfm2.CppMeshDynTri3D()
    dfm2.meshdyntri3d_initialize(dmesh, mesh.np_pos, mesh.np_elm)
    dmesh.check()

  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1.0])
    mesh = cad.mesh(0.2)
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

  def test2(self):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_2k.ply")
    dmesh = dfm2.CppMeshDynTri3D()
    dfm2.meshdyntri3d_initialize(dmesh, msh.np_pos, msh.np_elm)
    dmesh.check()
    for itr in range(1000):
      itri0 = random.randint(0, dmesh.ntri() - 1)
      iedge0 = random.randint(0, 2)
      dmesh.delete_tri_edge(itri0, iedge0)


class Test_MeshDynTri2D(unittest.TestCase):
  def test0(self):
    dmsh = dfm2.MeshDynTri2D()
    dmsh.meshing_loops([[-1,-1, +1,-1, +1,+1, -1,+1],
                        [-0.5, -0.5, +0.5, -0.5, 0.0, +0.5]],
                       edge_length=0.1)
