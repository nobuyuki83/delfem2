import sys, random
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw

def main_CppMeshDynTri2D_0():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = cad.mesh(0.1)
  dmesh = dfm2.CppMeshDynTri2D()
  dfm2.meshdyntri2d_initialize(dmesh,mesh.np_pos, mesh.np_elm)
  dmesh.check()
  dfm2.glfw.winDraw3d([dmesh,dfm2.AxisXYZ()],winsize=(400,300))


def main_CppMeshDynTri2D_1():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = cad.mesh(0.3)
  dmesh = dfm2.CppMeshDynTri2D()
  dfm2.meshdyntri2d_initialize(dmesh,mesh.np_pos, mesh.np_elm)
  dmesh.check()
  for itr in range(10):
    itri0 = random.randint(0,dmesh.ntri()-1)
    r0 = random.uniform(0.02, 0.98)
    r1 = random.uniform(0.01, 0.99-r0)
    ipo= dmesh.insert_point_elem(itri0,r0,r1)
    dmesh.delaunay_around_point(ipo)
    dmesh.check()
  dfm2.glfw.winDraw3d([dmesh],winsize=(400,300))


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
  dfm2.glfw.winDraw3d([dmesh],winsize=(400,300))


def main_MeshDynTri2D_0():
  dmsh = dfm2.MeshDynTri2D()
  dmsh.meshing_loops([[-1,-1, +1,-1, +1,+1, -1,+1],
                      [-0.5, -0.5, +0.5, -0.5, 0.0, +0.5]],
                     edge_length=0.2)
  dfm2.glfw.winDraw3d([dmsh])


def main_MeshDynTri2D_1():
  dmsh = dfm2.MeshDynTri2D()
  dmsh.meshing_loops([[-1,-1, +1,-1, +1,+1, -1,+1]],
                     edge_length=0.2)
  dfm2.glfw.winDraw3d([dmsh])
  res = dmsh.refine_EdgeLongerThan_InsideCircle(0.05, 0.0,0.0, 0.5)
  dfm2.glfw.winDraw3d([dmsh])


def main_MeshDynTri2D_2():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  dmsh = dfm2.MeshDynTri2D()
  dmsh.set_mesh(cad.mesh(0.1))
  fem = dfm2.FEM_Poisson(dmsh,source=1.0)
  npIdP = cad.points_edge([0,1,2,3], dmsh.np_pos)
  fem.value[npIdP] = 0.0
  fem.ls.bc[npIdP] = 1
  fem.solve()
  vis = dfm2.VisFEM_ColorContour(fem,name_color="value")
  dfm2.glfw.winDraw3d([vis,dmsh])
  #####
  res = dmsh.refine_EdgeLongerThan_InsideCircle(0.05, 0.0,0.0,0.5)
  fem.updated_mesh(res)
  dfm2.glfw.winDraw3d([vis,dmsh])
  #####
  fem.value[npIdP] = 0.0
  fem.ls.bc[npIdP] = 1
  fem.solve()
  dfm2.glfw.winDraw3d([vis,dmsh])


def main_MeshDynTri2D_3():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  dmsh = dfm2.MeshDynTri2D()
  dmsh.set_mesh(cad.mesh(0.1))
  fem = dfm2.FEM_Cloth(dmsh)
  npIdP = cad.points_edge([0], dmsh.np_pos)
  fem.ls.bc[npIdP,:] = 1
  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=dmsh.np_elm)
  dfm2.glfw.winDraw3d([fem,mesh2])
  #####
  mesh2 = None
  res = dmsh.refine_EdgeLongerThan_InsideCircle(0.05, 0.0,0.0,0.5)
  fem.updated_mesh(res)
  fem.ls.bc[npIdP,:] = 1
  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=dmsh.np_elm)
  dfm2.glfw.winDraw3d([fem,mesh2])
  '''  
  #####
  fem.vec_val[npIdP] = 0.0
  fem.ls.vec_bc[npIdP] = 1
  fem.solve()
  field = dfm2.Field(dmsh,val_color=fem.vec_val[:,0])
  dfm2.winDraw3d([field,dmsh])
  '''


if __name__ == "__main__":
  main_CppMeshDynTri2D_0()
  main_CppMeshDynTri2D_1()

  main_CppMeshDynTri3D()

  main_MeshDynTri2D_0()
  main_MeshDynTri2D_1()
  main_MeshDynTri2D_2()
  main_MeshDynTri2D_3()