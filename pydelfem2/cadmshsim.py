import numpy
import OpenGL.GL as gl

from .fem import FEM_Poisson, FEM_SolidLinearStatic, FEM_SolidLinearEigen, FEM_Diffuse
from .fem import PBD, PBD_Cloth
from .fem import FieldValueSetter
from .fem import VisFEM_ColorContour

from .cadmsh import CadMesh2D, cad_getPointsEdge, Mesh

from .gl.c_gl import setSomeLighting



class CadMesh2D_FEMPoisson(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_Poisson(source=1.0)
    self.vis = VisFEM_ColorContour(self.fem,name_color='value')
    self.vis.is_lighting = False
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
    npIdP = self.map_cad2msh.npIndEdge(self.list_cad_edge_fix)
    self.fem.value[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1
    self.fem.solve()


class CadMesh2D_FEMDiffuse(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_Diffuse()
    self.vis = VisFEM_ColorContour(self.fem,name_color='value')
    self.list_cad_edge_fix = [0,1,2,3]

  def draw(self):
    self.ccad.draw()
    self.vis.draw()

  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)

  def remesh(self):
    super().remesh()
    self.fem.updated_topology(self.dmsh)
    npIdP = cad_getPointsEdge(self.ccad, self.list_cad_edge_fix, self.dmsh.np_pos, 0.001)
    self.fem.value[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1


class CadMesh2D_FEMSolidLinearStatic(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_SolidLinearStatic()
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


class CadMesh2D_FEMSolidLinearEigen(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_SolidLinearEigen()
    self.vis = VisFEM_ColorContour(self.fem, name_disp="mode")
    self.remesh()

  def draw(self):
    self.ccad.draw()
    gl.glEnable(gl.GL_LIGHTING)
    setSomeLighting()
    gl.glColor3d(0.8,0.8,0.8)
    gl.glDisable(gl.GL_CULL_FACE)
    self.vis.draw()
#    self.msh25.draw()

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



class CadMesh2D_PBD(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.pbd = PBD()

  def draw(self):
    self.ccad.draw()
    gl.glDisable(gl.GL_LIGHTING)
    gl.glColor3d(0.8,0.8,0.8)
    self.vis_mesh.draw()

  def remesh(self):
    super().remesh()
    self.pbd.updated_topology(self.dmsh)
    npIdP = self.map_cad2msh.npIndEdge(0)
    self.pbd.vec_bc[npIdP] = 1
#    self.fvs = FieldValueSetter("0.3*sin(2*t)", self.pbd.vec_val, 0,
#                                mesh=self.dmsh, npIdP=npIdP, dt=self.pbd.dt)
    self.vis_mesh = Mesh(np_pos=self.pbd.vec_val,np_elm=self.dmsh.np_elm)

  def step_time(self):
    self.pbd.step_time()


class CadMesh2D_PBDCloth(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.pbd = PBD_Cloth()

  def draw(self):
    self.ccad.draw()
    gl.glDisable(gl.GL_LIGHTING)
    gl.glColor3d(0.8,0.8,0.8)
    self.vis_mesh.draw()

  def remesh(self):
    super().remesh()
    self.pbd.updated_topology(self.dmsh)
    npIdP = self.map_cad2msh.npIndEdge(3)
    self.pbd.bc[npIdP] = 1
    self.vis_mesh = Mesh(np_pos=self.pbd.vec_val,np_elm=self.dmsh.np_elm)

  def step_time(self):
    self.pbd.step_time()