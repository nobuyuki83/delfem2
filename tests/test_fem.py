import unittest, numpy, random
import PyDelFEM2 as dfm2

class Test_PBD(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    mesh,map_cad2mesh = cad.mesh(edge_len=0.2)
    pbd = dfm2.PBD()
    pbd.updated_topology(mesh)
    npIdP = cad.points_edge([0], mesh.np_pos)
    pbd.vec_bc[npIdP] = 1
    fvs = dfm2.FieldValueSetter("0.3*sin(20*t)", pbd.vec_val, 0,
                                mesh=mesh, npIdP=npIdP, dt=pbd.dt)
    for itr in range(100):
      fvs.step_time()
      pbd.step_time()

  def test_pbd_hex(self):
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
    for itr in range(100):
      fvs.step_time()
      pbd.step_time()


class Test_FEMPoission2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh,map_cad2msh = cad.mesh(0.02)
    fem = dfm2.FEM_Poisson(source=1.0)
    fem.updated_topology(msh)
    npIdP = cad.points_edge([0,1,2,3], msh.np_pos)
    fem.ls.bc[npIdP] = 1
    fem.solve()


class Test_FEMPoission3D(unittest.TestCase):
  def test1(self):
    sdf = dfm2.SDF()
    sdf.list_sdf.append(dfm2.CppSDF_Sphere(0.55, [-0.5, 0, 0], True))
    sdf.list_sdf.append(dfm2.CppSDF_Sphere(0.55, [+0.5, 0, 0], True))
    np_xyz, np_tet = dfm2.isosurface(sdf.list_sdf)
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


class Test_FEMDiffuse2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh,map_cad2msh = cad.mesh(0.02)
    fem = dfm2.FEM_Diffuse()
    fem.updated_topology(msh)
    npIdP = cad.points_edge([0, 1, 2, 3], msh.np_pos)
    fem.ls.bc[npIdP] = 1
    for itr in range(100):
      fem.step_time()


class Test_FemDiffuse3D(unittest.TestCase):
  def test1(self):
    sdf = dfm2.SDF()
    sdf.list_sdf.append(dfm2.CppSDF_Sphere(0.55, [-0.5, 0, 0], True))
    sdf.list_sdf.append(dfm2.CppSDF_Sphere(0.55, [+0.5, 0, 0], True))
    np_xyz, np_tet = dfm2.isosurface(sdf.list_sdf)
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
    for itr in range(100):
      fem.step_time()


class Test_FEMSolidLLinearStatic2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh,map_cad2msh = cad.mesh(0.02)
    fem = dfm2.FEM_SolidLinearStatic()
    fem.updated_topology(msh)
    npIdP = cad.points_edge([3], msh.np_pos)
    fem.ls.bc[npIdP,:] = 1
    fem.solve()


class Test_FEMSolidLLinearDynamic2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh,map_cad2msh = cad.mesh(0.02)
    fem = dfm2.FEM_SolidLinearDynamic()
    fem.param_gravity_x = -0.01
    fem.updated_topology(msh)
    npIdP = cad.points_edge([3], msh.np_pos)
    fem.ls.bc[npIdP, :] = 1
    for itr in range(100):
      fem.step_time()


class Test_FEMSorkes2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh,map_cad2msh = cad.mesh(0.02)
    fem = dfm2.FEM_StorksStatic2D(msh)
    npIdP0 = cad.points_edge([0,1,2,3], msh.np_pos)
    fem.ls.bc[npIdP0,0:2] = 1
    npIdP1 = cad.points_edge([2], msh.np_pos)
    fem.vec_val[npIdP1,0] = 1.0
    fem.solve()


class Test_FEMCloth(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh,map_cad2msh = cad.mesh(edge_len=0.05)
    fem = dfm2.FEM_Cloth()
    fem.updated_topology(msh)
    npIdP = cad.points_edge([2], msh.np_pos)
    fem.ls.bc[npIdP,0:3] = 1
    for itr in range(100):
      fem.step_time()




if __name__ == "__main__":
    unittest.main()