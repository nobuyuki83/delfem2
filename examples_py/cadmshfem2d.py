import sys
sys.path.append("../module_py")
import dfm2
import dfm2.glfw

def mesh():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  cmf = dfm2.CadMesh2D(cad,edge_length=0.1)
  dfm2.glfw.winDraw3d([cmf])

class CadMesh_Poisson(dfm2.CadMesh2D):

  def __init__(self,cad,edge_length:float):
    super().__init__(cad,edge_length)
    self.fem = dfm2.FEM_Poisson(self.msh,source=1.0)
    self.vis = dfm2.VisFEM_ColorContour(self.fem,name_color='value')
    self.remesh()

  def remesh(self):
    super().remesh()
    self.fem.updated_mesh()
    npIdP = self.cad.points_edge([0, 1, 2, 3], self.msh.np_pos)
    self.fem.value[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1
    self.fem.solve()
    
  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    self.fem.solve()

  def draw(self):
    self.cad.draw()
    self.vis.draw()

  def minmax_xyz(self):
    return self.msh.minmax_xyz()


def poisson():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  cmf = CadMesh_Poisson(cad,edge_length=0.1)
  dfm2.glfw.winDraw3d([cmf])

if __name__ == "__main__":
  mesh()
  poisson()