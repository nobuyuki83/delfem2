import numpy
import OpenGL.GL as gl


from .fem import FEM_Poisson, FEM_SolidLinearStatic, FEM_SolidLinearEigen
from .fem import VisFEM_ColorContour

from .cadmsh import CadMesh2D, cad_getPointsEdge, Mesh


class CadMesh2D_Poisson(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_Poisson(source=1.0)
    self.vis = VisFEM_ColorContour(self.fem,name_color='value')
    self.list_cad_edge_fix = [0,1,2,3]

  def draw(self):
    self.ccad.draw()
    self.vis.draw()

  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    if self.is_sync_mesh:
      self.fem.solve()

  def remesh(self):
    super().remesh()
    self.fem.updated_topology(self.dmsh)
    npIdP = cad_getPointsEdge(self.ccad, self.list_cad_edge_fix, self.dmsh.np_pos, 0.001)
    self.fem.value[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1
    self.fem.solve()


class CadMesh2D_SolidLinearStatic(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_SolidLinearStatic(gravity=[0, -0.1])
    self.vis = VisFEM_ColorContour(self.fem, name_disp="vec_val")
    self.list_cad_edge_fix = [3]
    self.remesh()

  def draw(self):
    self.ccad.draw()
    self.vis.draw()

  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    self.fem.solve()

  def remesh(self):
    super().remesh()
    self.fem.updated_topology(self.dmsh)
    npIdP = self.map_cad2msh.npIndEdge(3)
    self.fem.vec_val[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1
    self.fem.solve()


class CadMesh2D_SolidLinearEigen(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_SolidLinearEigen()
    self.vis = VisFEM_ColorContour(self.fem, name_disp="mode")
    self.remesh()

  def draw(self):
    self.ccad.draw()
    gl.glDisable(gl.GL_LIGHTING)
    gl.glColor3d(0.8,0.8,0.8)
    self.vis.draw()

  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    self.msh25.set_extrude(self.dmsh, 1)
    self.msh25.np_pos[:, 2] *= 0.05
    self.fem.solve()

  def remesh(self):
    super().remesh()
    self.msh25 = Mesh()
    self.msh25.set_extrude(self.dmsh, 1)
    self.msh25.np_pos[:, 2] *= 0.05
    self.fem.updated_topology(self.msh25)
    self.fem.ls.f[:] = numpy.random.uniform(-1, 1, self.msh25.np_pos.shape)
    self.fem.solve()