from .fem import FEM_Poisson, VisFEM_ColorContour
from .cadmsh import CadMesh2D, cad_getPointsEdge


class CadMesh2D_Poisson(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_Poisson(self.dmsh,source=1.0)
    self.vis = VisFEM_ColorContour(self.fem,name_color='value')
    self.remesh()

  def draw(self):
    self.ccad.draw()
    self.vis.draw()

  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    self.fem.solve()

  def remesh(self):
    super().remesh()
    self.fem.updated_mesh()
    npIdP = cad_getPointsEdge(self.ccad, [0, 1, 2, 3], self.dmsh.np_pos, 0.001)
    self.fem.value[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1
    self.fem.solve()