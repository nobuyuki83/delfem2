import numpy, os
import OpenGL.GL as gl
from .dfm2 import *

####################

class Mesh():
  def __init__(self,
               np_pos=numpy.ndarray((0,3),dtype=numpy.float64),
               np_elm=numpy.ndarray((0,3),dtype=numpy.int32),
               elem_type=None):
    print("PyMesh -- construct")
    assert type(np_pos) == numpy.ndarray
    assert type(np_elm) == numpy.ndarray
    assert np_pos.dtype == numpy.float64
    assert np_elm.dtype == numpy.int32
    self.color_face = [0.8, 0.8, 0.8, 1.0]
    self.is_draw_edge = True
    self.is_draw_face = True
    self.np_pos = np_pos
    self.np_elm = np_elm
    self.elem_type = elem_type

  def ndim(self) -> int:
    return self.np_pos.shape[1]

  def scale_xyz(self,scale:float):
    self.np_pos *= scale

  def minmax_xyz(self):
    if self.np_pos.shape[0] == 0:
      return [1,-1, 0,0, 0,0]
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
    if self.is_draw_face:
      gl.glColor4d(self.color_face[0], self.color_face[1], self.color_face[2], self.color_face[3])
      gl.glMaterialfv(gl.GL_FRONT_AND_BACK, gl.GL_DIFFUSE, self.color_face)
      draw_mesh_facenorm(self.np_pos,self.np_elm)
    if self.is_draw_edge:
      gl.glDisable(gl.GL_LIGHTING)
      gl.glLineWidth(1)
      gl.glColor3d(0,0,0)
      draw_mesh_edge(self.np_pos,self.np_elm)

  def subdiv(self):
    if self.elem_type == QUAD:
      np_pos1,np_quad1 = meshquad3d_subdiv(self.np_pos,self.np_elm)
      return Mesh(np_pos1,np_quad1,QUAD)
    if self.elem_type == HEX:
      np_pos1,np_quad1 = meshhex3d_subdiv(self.np_pos,self.np_elm)
      return Mesh(np_pos1,np_quad1,HEX)

  def meshtri2d(self,list_pos,list_elm):
    self.np_pos = numpy.array(list_pos, dtype=numpy.float32).reshape((-1, 2))
    self.np_elm = numpy.array(list_elm, dtype=numpy.int).reshape((-1, 3))

  def psup(self) -> numpy.ndarray:
    res = jarray_mesh_psup(self.np_elm, self.np_pos.shape[0])
    return res

  def write_obj(self,path_file) -> None:
    meshtri3d_write_obj(path_file, self.np_pos, self.np_elm)

  def read(self,path_file) -> None:
    if os.path.isfile(path_file):
      ext = path_file.rsplit(".", 1)[1]
      if ext == 'ply':
        self.np_pos, self.np_elm = meshtri3d_read_ply(path_file)
      if ext == 'obj':
        self.np_pos, self.np_elm = meshtri3d_read_obj(path_file)
      if ext == 'nas' or ext == 'bdf':
        self.np_pos, self.np_elm = meshtri3d_read_nastran(path_file)
      self.elem_type = TRI

#######################################

def mesh_grid(shape:list) -> Mesh:
  np_pos,np_quad = meshquad2d_grid(shape[0],shape[1])
  return Mesh(np_pos,np_quad,QUAD)

def mesh_cad(cad,len) -> Mesh:
  xy,tri = getMesh_cad(cad,len)
  mesh = Mesh(xy,tri,TRI)
  return mesh

###########################################################################

class MeshDynTri2D(Mesh):
  def __init__(self):
    super().__init__()
    self.elem_type = TRI
    self.np_pos = numpy.ndarray((0,2),dtype=numpy.float64)
    self.np_elm = numpy.ndarray((0,3),dtype=numpy.int32)
    self.dmsh = CppMeshDynTri2D()

  def set_mesh(self,msh=Mesh):
    assert msh.elem_type == TRI
    assert msh.np_pos.shape[1] == 2
    self.np_pos.resize(msh.np_pos.shape)
    self.np_elm.resize(msh.np_elm.shape)
    self.np_pos[:,:] = msh.np_pos
    self.np_elm[:,:] = msh.np_elm
    self.elem_type = msh.elem_type
    meshdyntri2d_initialize(self.dmsh, self.np_pos, self.np_elm)

  def meshing_loops(self,loops:list,edge_length:float):
    self.dmsh.meshing_loops(loops,edge_length)
    self.np_pos.resize((self.dmsh.npoint(),2))
    self.np_elm.resize((self.dmsh.ntri(),3))
    copyMeshDynTri2D(self.np_pos,self.np_elm, self.dmsh)

###########################################################################

class Grid3D:
  def __init__(self):
    self.vg = CppVoxelGrid()

  def add(self,ix,iy,iz):
    self.vg.add(ix,iy,iz)

  def mesh_quad3d(self) -> Mesh:
    list_xyz, list_tri = meshquad3d_voxelgrid(self.vg)
    np_pos = numpy.array(list_xyz, dtype=numpy.float64).reshape((-1, 3))
    np_elm = numpy.array(list_tri, dtype=numpy.int32).reshape((-1, 4))
    return Mesh(np_pos, np_elm, QUAD)

  def mesh_hex3d(self) -> Mesh:
    list_xyz, list_tri = meshhex3d_voxelgrid(self.vg)
    np_pos = numpy.array(list_xyz, dtype=numpy.float64).reshape((-1, 3))
    np_elm = numpy.array(list_tri, dtype=numpy.int32).reshape((-1, 8))
    return Mesh(np_pos, np_elm, HEX)



####################

class Cad2D():
  def __init__(self,list_xy=None):
    self.cad = CppCad2D()
    if not list_xy is None:
      self.cad.add_polygon(list_xy)

  def draw(self):
    self.cad.draw()

  def mouse(self,btn,action,mods,src,dir,view_height):
    self.cad.mouse(btn,action,mods,src,dir,view_height)

  def motion(self,src0,src1,dir):
    self.cad.motion(src0,src1,dir)

  def add_polygon(self,list_xy):
    self.cad.add_polygon(list_xy)

  def mesh(self,edge_len=0.05) -> Mesh:
    return mesh_cad(self.cad, edge_len)

  def points_edge(self, list_edge_index, np_xy, tolerance=0.01):
    return cad_getPointsEdge(self.cad,list_edge_index, np_xy, tolerance=tolerance)

  def getVertexXY_face(self,iface:int):
    list_xy_bound = self.cad.getVertexXY_face(0)
    return numpy.array(list_xy_bound).reshape([-1,2])

  def mvc(self,msh:Mesh):
    np_xy_bound = self.getVertexXY_face(0)
    W = mvc(msh.np_pos, np_xy_bound)
    assert W.shape[0] == msh.np_pos.shape[0]
    assert W.shape[1] == np_xy_bound.shape[0]
    return W


class CadMesh2D():
  def __init__(self,cad,edge_length:float):
    self.cad = cad
    self.cad.cad.is_draw_face = False
    self.edge_length = edge_length
    self.msh = cad.mesh(edge_len=self.edge_length)
    self.W = self.cad.mvc(self.msh)

  def draw(self):
    self.cad.draw()
    self.msh.draw()

  def mouse(self,btn,action,mods,src,dir,view_height):
    self.cad.mouse(btn,action,mods,src,dir,view_height)

  def motion(self,src0,src1,dir):
    self.cad.motion(src0,src1,dir)
    np_xy_bound = self.cad.getVertexXY_face(0)
    self.msh.np_pos = numpy.dot(self.W,np_xy_bound)
    max_asp,min_area = quality_meshTri2D(self.msh.np_pos,self.msh.np_elm)
    if max_asp > 5.0 or min_area < 0.0:
      self.remesh()

  def remesh(self):
    msh1 = self.cad.mesh(edge_len=self.edge_length)
    self.msh.np_pos = msh1.np_pos
    self.msh.np_elm = msh1.np_elm
    self.W = self.cad.mvc(self.msh)


#####################################################

class SDF():
  def __init__(self):
    self.list_sdf = []

  def draw(self):
    for sdf in self.list_sdf:
      sdf.draw()