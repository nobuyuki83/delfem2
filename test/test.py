import unittest, numpy, random
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


class Test_MathExpression(unittest.TestCase):
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

class Test_DynamicMesh(unittest.TestCase):
  def test0(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1.0])
    mesh = cad.mesh(0.1)
    dmesh = dfm2.CppMeshDynTri()
    dfm2.meshdyntri3d_initialize(dmesh, mesh.np_pos, mesh.np_elm)
    dmesh.check()

  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1.0])
    mesh = cad.mesh(0.2)
    dmesh = dfm2.CppMeshDynTri()
    dfm2.meshdyntri3d_initialize(dmesh, mesh.np_pos, mesh.np_elm)
    dmesh.check()
    for itr in range(100):
      itri0 = random.randint(0, dmesh.ntri() - 1)
      r0 = random.uniform(0.02, 0.98)
      r1 = random.uniform(0.01, 0.99 - r0)
      ipo = dmesh.insert_point_elem(itri0, r0, r1)
      dmesh.delaunay_around_point(ipo)
      dmesh.check()

  def main2(self):
    msh = dfm2.mesh_read("../test_inputs/bunny_2k.ply");
    dmesh = dfm2.CppMeshDynTri()
    dfm2.meshdyntri3d_initialize(dmesh, msh.np_pos, msh.np_elm)
    dmesh.check()
    for itr in range(1000):
      itri0 = random.randint(0, dmesh.ntri() - 1)
      iedge0 = random.randint(0, 2)
      dmesh.delete_tri_edge(itri0, iedge0)


class Test_PBD(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesh = cad.mesh(edge_len=0.2)
    pbd = dfm2.PBD(mesh)
    npIdP = cad.points_edge([0], mesh.np_pos)
    pbd.vec_bc[npIdP] = 1
    fvs = dfm2.FieldValueSetter("0.3*sin(20*t)", pbd.vec_val, 0,
                                mesh=mesh, npIdP=npIdP, dt=pbd.dt)
    for itr in range(100):
      fvs.step_time()
      pbd.step_time()

  def test_pbd_hex(voxelgrid):
    voxelgrid = dfm2.Grid3D()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(2, 0, 0)
    voxelgrid.add(1, 1, 0)
    msh = voxelgrid.mesh_hex3d()
    pbd = dfm2.PBD(msh)
    npIdP = numpy.array([0, 1, 2, 3], dtype=numpy.int32)
    pbd.vec_bc[npIdP] = 1
    fvs = dfm2.FieldValueSetter("0.4*sin(0.8*t)", pbd.vec_val, 1,
                                mesh=msh, npIdP=npIdP, dt=pbd.dt)
    for itr in range(100):
      fvs.step_time()
      pbd.step_time()


class Test_FEMPoission2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_Poisson(msh,source=1.0)
    npIdP = cad.points_edge([0,1,2,3], msh.np_pos)
    fem.ls.vec_bc[npIdP] = 1
    fem.solve()


class Test_FEMPoission3D(unittest.TestCase):
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


class Test_FEMSolidLLinearStatic2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_LinearSolidStatic(msh,gravity=[0,-0.1])
    npIdP = cad.points_edge([3], msh.np_pos)
    fem.ls.vec_bc[npIdP,:] = 1
    fem.solve()


class Test_FEMSolidLLinearDynamic2D(unittest.TestCase):
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