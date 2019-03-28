import numpy, os
import OpenGL.GL as gl
from .dfm2 import *

def normalize_rigmsh(rigmsh):
  aabb = rigmsh.aabb()
  center = aabb.center()
  rigmsh.translate(-center[0],-center[1],-center[2])
  rigmsh.scale(1.0/aabb.max_length())

####################

class Mesh():

  def __init__(self,
               np_pos=numpy.ndarray((0,3),dtype=numpy.float32),
               np_elm=numpy.ndarray((0,3),dtype=numpy.int)):
    print("PyMesh -- construct")
    assert type(np_pos) == numpy.ndarray
    assert type(np_elm) == numpy.ndarray
    self.color_face = [0.8, 0.8, 0.8, 1.0]
    self.is_draw_edge = True
    self.np_pos = np_pos
    self.np_elm = np_elm

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

  def psup(self):
    res = get_psup(self.np_elm, self.np_pos.shape[0])
    return res

#######################################

def mesh_read(path_file="") -> Mesh:
  if os.path.isfile(path_file):
    ext = path_file.rsplit(".", 1)[1]
    if ext == 'ply':
      list_xyz, list_tri = read_ply(path_file)
      np_pos = numpy.array(list_xyz, dtype=numpy.float32).reshape((-1, 3))
      np_elm = numpy.array(list_tri, dtype=numpy.int).reshape((-1, 3))
    if ext == 'obj':
      list_xyz, list_tri = read_obj(path_file)
      np_pos = numpy.array(list_xyz, dtype=numpy.float32).reshape((-1, 3))
      np_elm = numpy.array(list_tri, dtype=numpy.int).reshape((-1, 3))
#      print(self.np_pos.shape, self.np_elm.shape)
  return Mesh(np_pos,np_elm)

def mesh_voxelgrid(voxelgrid) -> Mesh:
  list_xyz, list_tri = getmesh_voxelgrid(voxelgrid)
  np_pos = numpy.array(list_xyz, dtype=numpy.float32).reshape((-1, 3))
  np_elm = numpy.array(list_tri, dtype=numpy.int).reshape((-1, 4))
  return Mesh(np_pos,np_elm)

def mesh_grid(shape) -> Mesh:
  h = shape[0]
  w = shape[1]
  list_xyz,list_elm = get_mesh_grid(h,w)
  msh0 = Mesh()
  msh0.np_pos = numpy.array(list_xyz, dtype=numpy.float32).reshape((-1, 2))
  msh0.np_elm = numpy.array(list_elm, dtype=numpy.int).reshape((-1, 4))
  return msh0

def mesh_cad(cad,len) -> Mesh:
  xy,tri = getMesh_cad(cad,len)
  mesh = Mesh(xy,tri)
  return mesh

#####################################################

class Field():
  def __init__(self,
               mesh: Mesh,
               val_color= None,
               val_disp=None):
    self.mesh = mesh
    self.val_color = val_color
    if type(val_color) == numpy.ndarray:
      self.draw_val_min = val_color.min()
      self.draw_val_max = val_color.max()
      self.color_mode = 'bcgyr'
    self.val_disp = val_disp

  def draw(self):
    if type(self.val_color) == numpy.ndarray:
      self.color_map = ColorMap(self.draw_val_min,self.draw_val_max,self.color_mode)
      drawField_colorMap(self.mesh.np_pos, self.mesh.np_elm,
                         self.val_color,
                         self.color_map)
    if type(self.val_disp) == numpy.ndarray:
      drawField_disp(self.mesh.np_pos, self.mesh.np_elm,
                     self.val_disp)

  def minmax_xyz(self):
    return self.mesh.minmax_xyz()

######################################################

class FEM_LinSys():
  def __init__(self,
               np:int, ndimval:int, pattern:tuple):
    # vectors
    self.vec_bc = numpy.zeros((np,ndimval), dtype=numpy.int32)
    self.vec_f = numpy.zeros((np,ndimval), dtype=numpy.float64)
    self.vec_x = numpy.zeros((np,ndimval), dtype=numpy.float64)

    # matrix
    self.mat = MatrixSquareSparse()
    self.mat.initialize(np, ndimval, True)
    psup_ind,psup = pattern[0],pattern[1]
    sortIndexedArray(psup_ind, psup)
    matrixSquareSparse_setPattern(self.mat, psup_ind, psup)

    # preconditioner
    self.mat_prec = PreconditionerILU()
    precond_ilu0(self.mat_prec, self.mat)

    self.conv_hist = []

  def SetZero(self):
    self.mat.setZero()
    self.vec_f[:,:] = 0.0

  def Solve(self):
    #### setting bc
    self.vec_f[self.vec_bc != 0] = 0.0
    matrixSquareSparse_setFixBC(self.mat, self.vec_bc)

    #### solving matrix
    self.mat_prec.set_value(self.mat)
    self.mat_prec.ilu_decomp()
    self.conv_hist = linsys_solve_pcg(self.vec_f, self.vec_x,
                                 0.0001, 100, self.mat, self.mat_prec)
    self.vec_x[self.vec_bc != 0] = 0.0

########################################


class FEM_LinearSolid2DStatic():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    self.ls = FEM_LinSys(np,2,pattern=mesh.psup())
    self.vec_val = numpy.zeros((np,2), dtype=numpy.float64)  # initial guess is zero

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_linearSolid2DStatic(self.ls.mat, self.ls.vec_f,
                                    1.0, 0.0, 1.0, 0.0, -0.1,
                                    self.mesh.np_pos, self.mesh.np_elm, self.vec_val)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x


class FEM_Poisson2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    self.ls = FEM_LinSys(np,1,mesh.psup())
    self.vec_val = numpy.zeros((np,1), dtype=numpy.float64)  # initial guess is zero

  def solve(self):
    self.ls.SetZero()

    mergeLinSys_poission2D(self.ls.mat, self.ls.vec_f,
                           1.0, 0.1,
                           self.mesh.np_pos, self.mesh.np_elm, self.vec_val)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x
