import numpy
from typing import List

from .c_core import TRI

from .fem import \
  FEM_ScalarPoisson, \
  FEM_ScalarDiffuse, \
  FEM_SolidLinearStatic, \
  FEM_SolidLinearEigen,\
  FEM_ShellPlateBendingMITC3_Eigen
from .fem import PBD, PBD_Cloth
from .fem import FieldValueSetter

from .cadmsh import CadMesh2D, cad_getPointsEdge, Mesh

from .gl.c_gl import setSomeLighting

import OpenGL.GL as gl
from .gl import VisFEM_ColorContour


class CadMesh2D_FEMPoisson(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_ScalarPoisson(source=1.0)
    self.vis = VisFEM_ColorContour(self.fem,name_color='value')
    self.vis.is_lighting = False
    self.list_cad_edge_fix = [0,1,2,3]

  def draw(self):
    self.ccad.draw()
    gl.glDisable(gl.GL_LIGHTING)
    self.vis.draw()

  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    if self.is_sync_mesh:
      self.fem.solve()

  def remesh(self):
    super().remesh()
    self.fem.updated_topology(self.dmsh)
    npIdP = self.mesher.points_on_edges(self.list_cad_edge_fix,self)
    self.fem.value[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1
    self.fem.solve()


class CadMesh2D_FEMDiffuse(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_ScalarDiffuse()
    self.vis = VisFEM_ColorContour(self.fem,name_color='value')
    self.vis.is_lighting = False
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
    self.vis.is_lighting = False
    self.list_cad_edge_fix = [3]

  def draw(self):
    self.ccad.draw()
    self.vis.draw()

  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    self.fem.solve()

  def remesh(self):
    super().remesh()
    self.fem.updated_topology(self.dmsh)
    npIdP = self.mesher.points_on_edges([3],self)
    self.fem.vec_val[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1
    self.fem.solve()


class CadMesh2D_FEMSolidLinearEigen(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_SolidLinearEigen()
    self.vis = VisFEM_ColorContour(self.fem, name_disp="mode")
    self.vis.is_lighting = True

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




class CadMesh2D_FEMShellPlateBendingMITC3Eigen(CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = FEM_ShellPlateBendingMITC3_Eigen()
    #self.vis = VisFEM_ColorContour(self.fem, name_disp="mode")
    #self.vis.is_lighting = True

  def draw(self):
    self.ccad.draw()
    gl.glEnable(gl.GL_LIGHTING)
    setSomeLighting()
    gl.glColor3d(0.8,0.8,0.8)
    gl.glDisable(gl.GL_CULL_FACE)
    self.vis_mesh.draw()
#    self.msh25.draw()

  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    self.fem.updated_geometry()
    self.fem.solve()
    # set visualization
    self.vis_mesh.np_pos[:, :2] = self.dmsh.np_pos
    self.vis_mesh.np_pos[:, 2] = self.fem.mode[:, 0]

  def remesh(self):
    """
    regenerate 2D mesh and data structure for linear system
    """
    super().remesh()
    self.fem.updated_topology(self.dmsh)
    self.fem.solve()
    self.vis_mesh = Mesh(numpy.zeros((self.dmsh.np_pos.shape[0], 3)),
                         self.dmsh.np_elm, TRI)
    self.vis_mesh.np_pos[:, :2] = self.dmsh.np_pos
    self.vis_mesh.np_pos[:, 2] = self.fem.mode[:, 0]

  def step_time(self):
    self.fem.solve()
    print(len(self.fem.ls.conv_hist), self.fem.freq_eigen)
    # set visualization
    self.vis_mesh.np_pos[:, :2] = self.dmsh.np_pos
    self.vis_mesh.np_pos[:, 2] = self.fem.mode[:, 0]



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
    npIdP = self.mesher.points_on_edges([0],self)
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
    self.listIndEdge_Fix = []
    self.list_seam = []
    self.pbd.param_gravity_x = 0
    self.pbd.param_gravity_y = 0
    self.pbd.param_gravity_z = 0

  def draw(self):
    self.ccad.draw()
    gl.glDisable(gl.GL_LIGHTING)
    gl.glColor3d(0.8,0.8,0.8)
    self.vis_mesh.draw()

  def remesh(self):
    super().remesh()
    self.pbd.updated_topology(self.dmsh)

    # fixbc
    npIdP = self.mesher.points_on_edges(self.listIndEdge_Fix,self)
    self.pbd.bc[npIdP] = 1

    # seam
    list_npIndPSeam = []
    for seam in self.list_seam:
      npIndP_Edge0a = self.mesher.points_on_one_edge(seam[0], True, self)
      npIndP_Edge0b = self.mesher.points_on_one_edge(seam[1], True, self)
      npIndP_Seam0 = numpy.vstack([npIndP_Edge0a, npIndP_Edge0b[::-1]]).transpose()
      list_npIndPSeam.append(npIndP_Seam0)
    if len(list_npIndPSeam) > 0:
      self.pbd.elems_seam = numpy.vstack(list_npIndPSeam).copy().astype(numpy.uint32)
    else:
      self.pbd.elems_seam = None

    # visualization
    self.vis_mesh = Mesh(np_pos=self.pbd.vec_val,np_elm=self.dmsh.np_elm)

  def step_time(self):
    self.pbd.step_time()