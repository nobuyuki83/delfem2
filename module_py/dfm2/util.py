import numpy, os
import OpenGL.GL as gl
from .dfm2 import *

def normalize_rigmsh(rigmsh):
  aabb = rigmsh.aabb()
  center = aabb.center()
  rigmsh.translate(-center[0],-center[1],-center[2])
  rigmsh.scale(1.0/aabb.max_length())


class Mesh():
  def __init__(self,path_file="",voxelgrid=None):
    self.color_face = [0.8, 0.8, 0.8, 1.0]
    self.is_draw_edge = True
    self.np_pos = numpy.empty((0,3),dtype=numpy.float32)
    self.np_elm = numpy.empty((0,3),dtype=numpy.int)
    if os.path.isfile(path_file):
      ext = path_file.rsplit(".",1)[1]
      if ext == 'ply':
        list_xyz, list_tri = read_ply(path_file)
        self.np_pos = numpy.array(list_xyz, dtype=numpy.float32).reshape((-1, 3))
        self.np_elm = numpy.array(list_tri, dtype=numpy.int).reshape((-1, 3))
      if ext == 'obj':
        list_xyz, list_tri = read_obj(path_file)
        self.np_pos = numpy.array(list_xyz, dtype=numpy.float32).reshape((-1, 3))
        self.np_elm = numpy.array(list_tri, dtype=numpy.int).reshape((-1, 3))
        print(self.np_pos.shape, self.np_elm.shape)
    if voxelgrid is not None:
      list_xyz, list_tri = getmesh_voxelgrid(voxelgrid)
      self.np_pos = numpy.array(list_xyz, dtype=numpy.float32).reshape((-1, 3))
      self.np_elm = numpy.array(list_tri, dtype=numpy.int).reshape((-1, 4))

  def ndim(self):
    return self.np_pos.shape[1]

  def scale_xyz(self,scale:float):
    self.np_pos *= scale

  def minmax_xyz(self):
    x_min = numpy.min(self.np_pos[:,0])
    x_max = numpy.max(self.np_pos[:,0])
    y_min = numpy.min(self.np_pos[:,1])
    y_max = numpy.max(self.np_pos[:,1])
    if self.np_pos.shape[1] >= 3:
      z_min = numpy.min(self.np_pos[:,2])
      z_max = numpy.max(self.np_pos[:,2])
    else:
      z_min = 0.0
      z_max = 0.0
    return [x_min,x_max, y_min,y_max, z_min,z_max]

  def draw(self):
    gl.glColor4d(self.color_face[0], self.color_face[1], self.color_face[2], self.color_face[3])
    gl.glMaterialfv(gl.GL_FRONT_AND_BACK, gl.GL_DIFFUSE, self.color_face)
    draw_mesh_facenorm(self.np_pos,self.np_elm)
    if self.is_draw_edge:
      gl.glDisable(gl.GL_LIGHTING)
      gl.glLineWidth(1)
      draw_mesh_edge(self.np_pos,self.np_elm)

  def subdiv(self):
    list_pos = numpy.ravel(self.np_pos).tolist()
    list_elm = numpy.ravel(self.np_elm).tolist()
    list_pos,list_elm = subdiv(list_pos,list_elm)
    msh0 = Mesh()
    msh0.np_pos = numpy.array(list_pos, dtype=numpy.float32).reshape((-1, 3))
    msh0.np_elm = numpy.array(list_elm, dtype=numpy.int).reshape((-1, 4))
    return msh0

  def meshtri2d(self,list_pos,list_elm):
    self.np_pos = numpy.array(list_pos, dtype=numpy.float32).reshape((-1, 2))
    self.np_elm = numpy.array(list_elm, dtype=numpy.int).reshape((-1, 3))

def grid_mesh(shape) -> Mesh:
  h = shape[0]
  w = shape[1]
  list_xyz,list_elm = get_mesh_grid(h,w)
  msh0 = Mesh()
  msh0.np_pos = numpy.array(list_xyz, dtype=numpy.float32).reshape((-1, 2))
  msh0.np_elm = numpy.array(list_elm, dtype=numpy.int).reshape((-1, 4))
  return msh0


class Field():
  def __init__(self,
               mesh: Mesh,
               val: numpy.ndarray):
    self.mesh = mesh
    self.val = val
    self.draw_val_min = val.min()
    self.draw_val_max = val.max()
    self.color_mode = 'bcgyr'


  def draw(self):
    self.color_map = ColorMap(self.draw_val_min,self.draw_val_max,self.color_mode)
    draw_field(self.mesh.np_pos, self.mesh.np_elm,
               self.val,
               self.color_map)

