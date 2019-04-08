import unittest
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


class Test_FEMPoission2D(unittest.TestCase):
  def test1(self):
    cad = dfm2.Cad2D(list_xy=[-1, -1, +1, -1, +1, +1, -1, +1])
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_Poisson(msh,source=1.0)
    npIdP = cad.points_edge([0,1,2,3], msh.np_pos)
    fem.ls.vec_bc[npIdP] = 1
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
    msh = cad.mesh(0.02)
    fem = dfm2.FEM_Cloth(msh)
    npIdP = cad.points_edge([2], msh.np_pos)
    fem.ls.vec_bc[npIdP,0:3] = 1
    for itr in range(100):
      fem.step_time()

if __name__ == "__main__":
    unittest.main()