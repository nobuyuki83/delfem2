import sys, random
sys.path.append("../module_py")
import dfm2

def main0():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = cad.mesh(0.1)
  dmesh = dfm2.CppMeshDynTri()
  dfm2.meshdyntri3d_initialize(dmesh,mesh.np_pos, mesh.np_elm)
  dmesh.check()
  dfm2.winDraw3d([dmesh],winsize=(400,300))

def main1():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = cad.mesh(0.2)
  dmesh = dfm2.CppMeshDynTri()
  dfm2.meshdyntri3d_initialize(dmesh,mesh.np_pos, mesh.np_elm)
  dmesh.check()
  for itr in range(100):
    itri0 = random.randint(0,dmesh.ntri()-1)
    r0 = random.uniform(0.02, 0.98)
    r1 = random.uniform(0.01, 0.99-r0)
    ipo= dmesh.insert_point_elem(itri0,r0,r1)
    dmesh.delaunay_around_point(ipo)
    dmesh.check()
  dfm2.winDraw3d([dmesh],winsize=(400,300))

def main2():
  msh = dfm2.mesh_read("../test_inputs/bunny_2k.ply");
  dmesh = dfm2.CppMeshDynTri()
  dfm2.meshdyntri3d_initialize(dmesh,msh.np_pos, msh.np_elm)
  dmesh.check()
  for itr in range(1000):
    itri0 = random.randint(0,dmesh.ntri()-1)
    iedge0 = random.randint(0,2)
    dmesh.delete_tri_edge(itri0,iedge0)
  dfm2.winDraw3d([dmesh],winsize=(400,300))

if __name__ == "__main__":
  main0()
  main1()
  main2()
