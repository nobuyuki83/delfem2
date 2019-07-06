import numpy

from ..cadmsh import Mesh

from ..c_core import CppMeshDynTri2D, CppMeshDynTri3D, CppCad2D, CppSDF_Sphere
from ..c_core import TRI, QUAD, LINE

from .c_gl import CppGPUSampler, depth_buffer
from .c_gl import ColorMap
from .c_gl import get_texture, setSomeLighting
from .c_gl import cppDrawEdge_CppMeshDynTri2D, cppDrawEdge_CppMeshDynTri3D, cppDraw_CppCad2D
from .c_gl import draw_mesh_facenorm, draw_mesh_edge
from .c_gl import cppDrawSphere

from ._gl import AxisXYZ, Camera
from ._gl import getOpenglInfo, screenUnProjection, screenUnProjectionDirection

import OpenGL.GL as gl


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
    assert msh.elem_type == TRI or msh.elem_type == QUAD or LINE
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

def draw_CppSDF_Sphere(self):
  c = self.cent
  gl.glDisable(gl.GL_LIGHTING)
  gl.glColor3d(1,0,0,)
  cppDrawSphere(32,32,self.rad,c[0],c[1],c[2])


setattr(Mesh,"draw",draw_Mesh)
setattr(Mesh,"is_draw_edge",True)
setattr(Mesh,"is_draw_face",True)
setattr(Mesh,"color_face", [0.8, 0.8, 0.8, 1.0])

setattr(CppMeshDynTri2D,"draw",draw_CppMeshDynTri2D)
setattr(CppMeshDynTri3D,"draw",draw_CppMeshDynTri3D)
setattr(CppCad2D,       "draw",draw_CppCad2D)
setattr(CppSDF_Sphere,  "draw",draw_CppSDF_Sphere)