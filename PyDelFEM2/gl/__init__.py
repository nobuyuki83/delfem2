####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import OpenGL.GL as gl
import numpy

from ..cadmsh import Mesh

from ..c_core import CppMeshDynTri2D, CppMeshDynTri3D, CppCad2D
from ..c_core import TRI, QUAD, LINE
from ..c_core import CppSDF3, CppSDF3_Sphere
#from ..c_core import RigidBodyAssembly_Static
from ..c_core import \
  write_vtk_meshpoint, \
  write_vtk_meshelem, \
  write_vtk_pointscalar, \
  write_vtk_pointvector

from .c_gl import CppGPUSampler, depth_buffer, color_buffer_4float
from .c_gl import ColorMap
from .c_gl import get_texture, setSomeLighting, setup_glsl
from .c_gl import cppDrawEdge_CppMeshDynTri2D, cppDrawEdge_CppMeshDynTri3D, cppDraw_CppCad2D
from .c_gl import draw_mesh_facenorm, draw_mesh_edge
from .c_gl import cppDrawSphere
from .c_gl import drawField_colorMap, drawField_disp, drawField_hedgehog

from ._gl import AxisXYZ, Camera, CAMERA_ROT_MODE
from ._gl import getOpenglInfo, screenUnProjection, screenUnProjectionDirection


class VisFEM_Hedgehog():
  def __init__(self, fem,
               name_vector=""):
    self.fem = fem
    self.name_vector = name_vector

  def draw(self):
    mesh = self.fem.mesh
    if hasattr(self.fem, self.name_vector):
      npVector = getattr(self.fem, self.name_vector)
      assert type(npVector) == numpy.ndarray
      assert npVector.ndim  == 2
      ndim = self.fem.mesh.np_pos.shape[1]
      gl.glDisable(gl.GL_LIGHTING)
      gl.glColor3d(0, 0, 0)
      drawField_hedgehog(self.fem.mesh.np_pos, npVector[:,0:ndim], 1.0)

class VisFEM_ColorContour():
  def __init__(self, fem,
               name_color="",
               name_disp="",
               idim = 0):
    self.fem = fem
    ####
    self.name_color = name_color
    self.idim = idim
    self.color_mode = 'bcgyr'
    self.color_min = 0.0
    self.color_max = 0.3
    ####
    self.name_disp = name_disp
    self.disp_mode = 'disp'
    self.is_lighting = True

  def minmax_xyz(self):
    return self.fem.mesh.minmax_xyz()

  def set_color_minmax(self):
    if hasattr(self.fem, self.name_color):
      npColor = getattr(self.fem, self.name_color)
      assert type(npColor) == numpy.ndarray
      assert npColor.ndim == 2 and self.idim < npColor.shape[1]
      self.color_min = npColor[:,self.idim].min()
      self.color_max = npColor[:,self.idim].max()

  def draw(self):
    if not hasattr(self.fem, 'mesh'):
      return
    mesh = self.fem.mesh
    if self.is_lighting:
      gl.glEnable(gl.GL_LIGHTING)

    if hasattr(self.fem, self.name_color):
      npColor = getattr(self.fem, self.name_color)
      assert type(npColor) == numpy.ndarray
      assert npColor.ndim == 2 and self.idim < npColor.shape[1]
      self.color_map = ColorMap(self.color_min,self.color_max,self.color_mode)
      drawField_colorMap(mesh.np_pos, mesh.np_elm,
                         npColor[:,self.idim],
                         self.color_map)

    if hasattr(self.fem, self.name_disp):
      npDisp = getattr(self.fem, self.name_disp)
      assert type(npDisp) == numpy.ndarray
      if self.disp_mode == 'disp':
        drawField_disp(mesh.np_pos, mesh.np_elm, mesh.elem_type,
                       npDisp)
      if self.disp_mode == 'hedgehog':
        gl.glDisable(gl.GL_LIGHTING)
        gl.glColor3d(0,0,0)
        drawField_hedgehog(mesh.np_pos, self.val_disp, 1.0)

  def write_vtk(self, path_vtk, message=""):
    mesh = self.fem.mesh
    write_vtk_meshpoint(path_vtk,"foobar", mesh.np_pos)
    write_vtk_meshelem(path_vtk, mesh.np_elm, mesh.elem_type)
    open(path_vtk, "a+").write("POINT_DATA {0}\n".format(mesh.np_pos.shape[0]))
    if hasattr(self.fem, self.name_color):
      npColor = getattr(self.fem, self.name_color)
      assert type(npColor) == numpy.ndarray
      if npColor.ndim == 1 or (npColor.ndim == 2 and npColor.shape[1] == 1):
        write_vtk_pointscalar(path_vtk, npColor)
    elif hasattr(self.fem, self.name_disp):
      npDisp = getattr(self.fem, self.name_disp)
      assert type(npDisp) == numpy.ndarray
      if npDisp.ndim == 2 and npDisp.shape[1] == mesh.np_pos.shape[1]:
        write_vtk_pointvector(path_vtk, npDisp)



class GLBufferMesh():
  def __init__(self,mesh=None,is_normal=True):
    self.vbo_pos = 0
    self.vbo_nrm = 0
    self.ebo_elm = 0
    self.ndim = 0
    self.size_elem = 0
    self.gl_elem_type = gl.GL_NONE

    if isinstance(mesh,Mesh):
      self.set_mesh(mesh,is_normal)

  def release_buffer(self):
    if gl.glIsBuffer(self.vbo_pos):
      gl.glDeleteBuffers(1, [self.vbo_pos])
    self.vbo_pos = 0
    ##
    if gl.glIsBuffer(self.vbo_nrm):
      gl.glDeleteBuffers(1, [self.vbo_nrm])
    self.vbo_nrm = 0
    ##
    if gl.glIsBuffer(self.ebo_elm):
      gl.glDeleteBuffers(1, [self.ebo_elm])
    self.ebo_elm = 0

  def set_mesh(self,msh:Mesh, is_normal:bool):
    assert msh.elem_type == TRI or msh.elem_type == QUAD or  msh.elem_type ==LINE
    ####
    self.release_buffer()
    ####
    self.vbo_pos = vbo_array(msh.np_pos)
    self.ebo_elm = ebo_array(msh.np_elm)
    self.ndim = msh.np_pos.shape[1]
    self.size_elem = msh.np_elm.size
    if msh.elem_type == TRI:
      self.gl_elem_type = gl.GL_TRIANGLES
    elif msh.elem_type == LINE:
      self.gl_elem_type = gl.GL_LINES
    else:
      assert(0)

    if is_normal and (msh.elem_type == TRI or msh.elem_type == QUAD):
      nrm = msh.normal_vtx()
      self.vbo_nrm = vbo_array(nrm)

  def draw(self):
    assert gl.glIsBuffer(self.vbo_pos)
    assert gl.glIsBuffer(self.ebo_elm)
    ####
    gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self.vbo_pos)
    gl.glVertexPointer(self.ndim, gl.GL_DOUBLE, 0, None)
    if gl.glIsBuffer(self.vbo_nrm):
      gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
      gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self.vbo_nrm)
      gl.glNormalPointer(gl.GL_DOUBLE, 0, None)
    gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, self.ebo_elm)
    gl.glDrawElements(self.gl_elem_type, self.size_elem, gl.GL_UNSIGNED_INT, None)
    gl.glDisableClientState(gl.GL_VERTEX_ARRAY)


def draw_Mesh(self):
  if self.is_draw_face:
    gl.glColor4d(self.color_face[0], self.color_face[1], self.color_face[2], self.color_face[3])
    gl.glMaterialfv(gl.GL_FRONT_AND_BACK, gl.GL_DIFFUSE, self.color_face)
    draw_mesh_facenorm(self.np_pos,self.np_elm, self.elem_type)

  if self.is_draw_edge:
    gl.glDisable(gl.GL_LIGHTING)
    gl.glLineWidth(1)
    gl.glColor3d(0,0,0)
    draw_mesh_edge(self.np_pos,self.np_elm, self.elem_type)


def vbo_array(aXY:numpy.ndarray) -> int:
  assert aXY.dtype == numpy.float64
  vbo = gl.glGenBuffers(1)
  gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)
  gl.glBufferData(gl.GL_ARRAY_BUFFER,
                  aXY.size * 8,
                  (gl.ctypes.c_double * aXY.size)(*aXY.flatten()),
                  gl.GL_STATIC_DRAW)
  gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
  return vbo


def ebo_array(aElm:numpy.ndarray) -> int:
  ebo = gl.glGenBuffers(1)
  gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, ebo)
  gl.glBufferData(gl.GL_ELEMENT_ARRAY_BUFFER,
                  aElm.size*4,
                  (gl.ctypes.c_uint32 * aElm.size)(*aElm.flatten()),
                  gl.GL_STATIC_DRAW)
  gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, 0)
  return ebo


def draw_CppMeshDynTri2D(self:CppMeshDynTri2D):
  cppDrawEdge_CppMeshDynTri2D(self)

def draw_CppMeshDynTri3D(self:CppMeshDynTri3D):
  cppDrawEdge_CppMeshDynTri3D(self)

def draw_CppCad2D(self:CppCad2D):
  cppDraw_CppCad2D(self)

def draw_CppSDF3_Sphere(self):
  c = self.cent
  gl.glDisable(gl.GL_LIGHTING)
  gl.glColor3d(1,0,0,)
  cppDrawSphere(32,32,self.rad,c[0],c[1],c[2])


#def draw_CppRigidBodyAssembly_Static(self:RigidBodyAssembly_Static):
#  drawRigidBodyAssemblyStatic(self)



setattr(Mesh,"draw",draw_Mesh)
setattr(Mesh,"is_draw_edge",True)
setattr(Mesh,"is_draw_face",True)
setattr(Mesh,"color_face", [0.8, 0.8, 0.8, 1.0])

setattr(CppMeshDynTri2D,"draw",draw_CppMeshDynTri2D)
setattr(CppMeshDynTri3D,"draw",draw_CppMeshDynTri3D)
setattr(CppCad2D,       "draw",draw_CppCad2D)
setattr(CppSDF3_Sphere, "draw",draw_CppSDF3_Sphere)
#setattr(RigidBodyAssembly_Static, "draw", draw_CppRigidBodyAssembly_Static)