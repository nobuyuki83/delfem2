import unittest, numpy
import sys
sys.path.append("../module_py")
import dfm2

class Test_Cad2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])

  def test2(self):
    cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])

  def test3(self):
    cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
    msh = cad.mesh(0.02)
    self.assertEqual(msh.np_pos.shape[1],2)
    W = cad.mvc(msh)
    self.assertEqual(W.ndim,2)
    self.assertEqual(W.shape[0],msh.np_pos.shape[0])
    self.assertEqual(W.shape[1],cad.getVertexXY_face(0).shape[0])
    self.assertLess(numpy.linalg.norm(W.sum(axis=1)-numpy.ones((W.shape[0]))),1.0e-3)


class TetMathExpression(unittest.TestCase):
  def test1(self):
    mee = dfm2.MathExpressionEvaluator()
    mee.set_key("x",3)
    mee.set_expression("x+3")
    mee.set_key("x",6)
    self.assertLess( mee.eval()-9.0, 1.0e-30 )

class Test_Mesh(unittest.TestCase):
  def test1(self):
    msh = dfm2.mesh_read("../test_inputs/bunny_2k.ply")
    self.assertIsNot(msh,None)

  def test2(self):
    msh = dfm2.mesh_read("../test_inputs/bunny_1k.obj");
    self.assertIsNot(msh,None)

  def test3(self):
    msh = dfm2.mesh_grid((32,64))
    self.assertIsNot(msh,None)

  def test4(self):
    voxelgrid = dfm2.VoxelGrid()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = dfm2.mesh_voxelgrid(voxelgrid=voxelgrid)
    msh = msh.subdiv()
    msh = msh.subdiv()
    msh = msh.subdiv()
    self.assertIsNot(msh,None)


class Test_FEMPoission2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_Poisson(msh,source=1.0)
    npIdP = cad.points_edge([0,1,2,3], msh.np_pos)
    fem.ls.vec_bc[npIdP] = 1
    fem.solve()


class TetFEMPoission3D(unittest.TestCase):
  def test1(self):
    sdf = dfm2.SDF()
    sdf.list_sdf.append(dfm2.SDF_Sphere(0.55, [-0.5, 0, 0], True))
    sdf.list_sdf.append(dfm2.SDF_Sphere(0.55, [+0.5, 0, 0], True))
    np_xyz, np_tet = dfm2.isosurface(sdf.list_sdf)
    msh = dfm2.Mesh(np_xyz, np_tet, dfm2.TET)
    npIdP0 = numpy.where(msh.np_pos[:,0]>+1)
    npIdP1 = numpy.where(msh.np_pos[:,0]<-1)
    fem = dfm2.FEM_Poisson(msh)
    fem.ls.vec_bc[npIdP0] = 1
    fem.ls.vec_bc[npIdP1] = 2
    fem.vec_val[:] = 0.5
    fem.vec_val[npIdP0] = 0.0
    fem.vec_val[npIdP1] = 1.0
    fem.solve()


class Test_FEMDiffuse2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_Diffuse(msh, source=1.0)
    npIdP = cad.points_edge([0, 1, 2, 3], msh.np_pos)
    fem.ls.vec_bc[npIdP] = 1
    for itr in range(100):
      fem.step_time()


class Test_FemDiffuse3D(unittest.TestCase):
  def test1(self):
    sdf = dfm2.SDF()
    sdf.list_sdf.append(dfm2.SDF_Sphere(0.55, [-0.5, 0, 0], True))
    sdf.list_sdf.append(dfm2.SDF_Sphere(0.55, [+0.5, 0, 0], True))
    np_xyz, np_tet = dfm2.isosurface(sdf.list_sdf)
    msh = dfm2.Mesh(np_xyz, np_tet, dfm2.TET)
    npIdP0 = numpy.where(msh.np_pos[:,0]>+1)
    npIdP1 = numpy.where(msh.np_pos[:,0]<-1)
    fem = dfm2.FEM_Diffuse(msh)
    fem.ls.vec_bc[npIdP0] = 1
    fem.ls.vec_bc[npIdP1] = 2
    fem.vec_val[:] = 0.5
    fem.vec_val[npIdP0] = 0.0
    fem.vec_val[npIdP1] = 1.0
    for itr in range(100):
      fem.step_time()


class TestFEM_SolidLLinearStatic2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_LinearSolidStatic(msh,gravity=[0,-0.1])
    npIdP = cad.points_edge([3], msh.np_pos)
    fem.ls.vec_bc[npIdP,:] = 1
    fem.solve()


class TestFEM_SolidLLinearDynamic2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_LinearSolidDynamic(msh, gravity=[0, -0.1])
    npIdP = cad.points_edge([3], msh.np_pos)
    fem.ls.vec_bc[npIdP, :] = 1
    for itr in range(100):
      fem.step_time()


class Test_FEMSorkes2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_StorksStatic2D(msh)
    npIdP0 = cad.points_edge([0,1,2,3], msh.np_pos)
    fem.ls.vec_bc[npIdP0,0:2] = 1
    npIdP1 = cad.points_edge([2], msh.np_pos)
    fem.vec_val[npIdP1,0] = 1.0
    fem.solve()


class Test_FEMCloth(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(edge_len=0.05)
    fem = dfm2.FEM_Cloth(msh)
    npIdP = cad.points_edge([2], msh.np_pos)
    fem.ls.vec_bc[npIdP,0:3] = 1
    for itr in range(100):
      fem.step_time()


if __name__ == "__main__":
    unittest.main()