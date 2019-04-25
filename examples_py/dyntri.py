import sys, random
sys.path.append("../module_py")
import dfm2

def main_CppMeshDynTri2D_0():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = cad.mesh(0.1)
  dmesh = dfm2.CppMeshDynTri2D()
  dfm2.meshdyntri2d_initialize(dmesh,mesh.np_pos, mesh.np_elm)
  dmesh.check()
  dfm2.winDraw3d([dmesh,dfm2.AxisXYZ()],winsize=(400,300))


def main_CppMeshDynTri2D_1():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = cad.mesh(0.3)
  dmesh = dfm2.CppMeshDynTri2D()
  dfm2.meshdyntri2d_initialize(dmesh,mesh.np_pos, mesh.np_elm)
  dmesh.check()
  for itr in range(1000):
    itri0 = random.randint(0,dmesh.ntri()-1)
    r0 = random.uniform(0.02, 0.98)
    r1 = random.uniform(0.01, 0.99-r0)
    ipo= dmesh.insert_point_elem(itri0,r0,r1)
    dmesh.delaunay_around_point(ipo)
    dmesh.check()
  dfm2.winDraw3d([dmesh],winsize=(400,300))


def main_MeshDynTri2D_0():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  dmsh = dfm2.MeshDynTri2D(cad.mesh(0.1))
  fem = dfm2.FEM_Poisson(dmsh,source=1.0)
  npIdP = cad.points_edge([0,1,2,3], dmsh.np_pos)
  fem.ls.vec_bc[npIdP] = 1
  fem.solve()
  field = dfm2.Field(dmsh,val_color=fem.vec_val[:,0])
  dfm2.winDraw3d([field,dmsh])
  ####



def main_CppMeshDynTri3D():
  msh = dfm2.Mesh()
  msh.read("../test_inputs/bunny_2k.ply")
  dmesh = dfm2.CppMeshDynTri3D()
  dfm2.meshdyntri3d_initialize(dmesh, msh.np_pos, msh.np_elm)
  dmesh.check()
  for itr in range(500):
    itri0 = random.randint(0,dmesh.ntri()-1)
    iedge0 = random.randint(0,2)
    dmesh.delete_tri_edge(itri0,iedge0)
  dfm2.winDraw3d([dmesh],winsize=(400,300))

if __name__ == "__main__":
  main_CppMeshDynTri2D_0()
  main_CppMeshDynTri2D_1()

  main_CppMeshDynTri3D()

  main_MeshDynTri2D_0()
