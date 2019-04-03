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
               np_pos=numpy.ndarray((0,3),dtype=numpy.float64),
               np_elm=numpy.ndarray((0,3),dtype=numpy.int32)):
    print("PyMesh -- construct")
    assert type(np_pos) == numpy.ndarray
    assert type(np_elm) == numpy.ndarray
    assert np_pos.dtype == numpy.float64
    assert np_elm.dtype == numpy.int32
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
    res = psup_mesh(self.np_elm, self.np_pos.shape[0])
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
               val_disp=None,
               disp_mode = 'disp'):
    self.mesh = mesh
    self.val_color = val_color
    self.disp_mode = disp_mode
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
      if self.disp_mode == 'disp':
        drawField_disp(self.mesh.np_pos, self.mesh.np_elm,
                       self.val_disp)
      if self.disp_mode == 'hedgehog':
        gl.glDisable(gl.GL_LIGHTING)
        gl.glColor3d(0,0,0)
        drawField_hedgehog(self.mesh.np_pos, self.val_disp, 1.0)

  def minmax_xyz(self):
    return self.mesh.minmax_xyz()

######################################################

class FEM_LinSys():
  def __init__(self,
               np:int, ndimval:int, pattern:tuple):
    self.nitr = 1000
    self.conv_ratio = 1.0e-4

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

  def Solve(self,is_asymmetric=False):
    #### setting bc
    self.vec_f[self.vec_bc != 0] = 0.0
    matrixSquareSparse_setFixBC(self.mat, self.vec_bc)

    #### solving matrix
    self.mat_prec.set_value(self.mat)
    self.mat_prec.ilu_decomp()
    if not is_asymmetric:
      self.conv_hist = linsys_solve_pcg(self.vec_f, self.vec_x,
                                        self.conv_ratio, self.nitr,
                                        self.mat, self.mat_prec)
    else:
      self.conv_hist = linsys_solve_bicgstab(self.vec_f, self.vec_x,
                                             self.conv_ratio, self.nitr,
                                             self.mat, self.mat_prec)
    self.vec_x[self.vec_bc != 0] = 0.0

########################################


class FEM_LinearSolidStatic2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    ndimval = 2
    np = mesh.np_pos.shape[0]
    self.ls = FEM_LinSys(np,ndimval,pattern=mesh.psup())
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_linearSolidStatic2D(self.ls.mat, self.ls.vec_f,
                                    1.0, 0.0, 1.0, 0.0, -0.1,
                                    self.mesh.np_pos, self.mesh.np_elm, self.vec_val)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x


class FEM_LinearSolidDynamic2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    ndimval = 2
    self.ls = FEM_LinSys(np,ndimval,pattern=mesh.psup())
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_acc = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.dt = 0.1
    self.gamma_newmark = 0.6
    self.beta_newmark = 0.36

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_linearSolidDynamic2D(self.ls.mat, self.ls.vec_f,
                                     1.0, 0.0, 1.0, 0.0, -0.1,
                                     self.dt, self.gamma_newmark, self.beta_newmark,
                                     self.mesh.np_pos, self.mesh.np_elm,
                                     self.vec_val, self.vec_velo, self.vec_acc)
    self.ls.Solve()
    self.vec_val += (self.dt)*self.vec_velo + (0.5*self.dt*self.dt)*self.vec_acc + (self.dt*self.dt*self.beta_newmark)*self.ls.vec_x
    self.vec_velo += (self.dt*self.gamma_newmark)*self.ls.vec_x + (self.dt)*self.vec_acc
    self.vec_acc += self.ls.vec_x

  def step_time(self):
    self.solve()


class FEM_Cloth():
  def __init__(self,mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    ndimval = 3
    ####
    self.np_quad = elemQuad_dihedralTri(mesh.np_elm, np)
    ####
    self.ls = FEM_LinSys(np,ndimval,pattern=psup_mesh(self.np_quad, np))
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.dt = 0.1
    self.vec_val[:,:2] = mesh.np_pos
    self.vec_val[:,2] = 1.0

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_cloth(self.ls.mat, self.ls.vec_f,
                      10.0, 500.0, self.dt,
                      self.mesh.np_pos, self.mesh.np_elm,
                      self.np_quad,
                      self.vec_val)
    mergeLinSys_massPoint(self.ls.mat, self.ls.vec_f,
                          1.0, self.dt,
                          [0,0,-1],
                          self.vec_val, self.vec_velo)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x
    self.vec_velo = (1.0/self.dt)*self.ls.vec_x

  def step_time(self):
    self.solve()


class FEM_Poisson2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    ndimval = 1
    self.ls = FEM_LinSys(np,ndimval,mesh.psup())
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_poission2D(self.ls.mat, self.ls.vec_f,
                           1.0, 0.1,
                           self.mesh.np_pos, self.mesh.np_elm,
                           self.vec_val)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x


class FEM_Diffuse2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    ndimval = 1
    self.ls = FEM_LinSys(np,ndimval,mesh.psup())
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.dt = 0.01
    self.gamma_newmark = 0.6

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_diffuse2D(self.ls.mat, self.ls.vec_f,
                          1.0, 1.0, 1.0,
                          self.dt, self.gamma_newmark,
                          self.mesh.np_pos, self.mesh.np_elm,
                          self.vec_val, self.vec_velo)
    self.ls.Solve()
    self.vec_val += (self.ls.vec_x)*(self.dt*self.gamma_newmark) + (self.vec_velo)*self.dt
    self.vec_velo += self.ls.vec_x

  def step_time(self):
    self.solve()


class FEM_StorksStatic2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    ndimval = 3
    self.ls = FEM_LinSys(np,ndimval,mesh.psup())
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_storksStatic2D(self.ls.mat, self.ls.vec_f,
                               1.0, 0.0, 0.0,
                               self.mesh.np_pos, self.mesh.np_elm,
                               self.vec_val)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x

  def step_time(self):
    self.solve()


class FEM_StorksDynamic2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    ndimval = 3
    self.ls = FEM_LinSys(np,ndimval,mesh.psup())
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.dt = 0.005
    self.gamma_newmark = 0.6

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_storksDynamic2D(self.ls.mat, self.ls.vec_f,
                                1.0, 1.0, 0.0, 0.0,
                                self.dt, self.gamma_newmark,
                                self.mesh.np_pos, self.mesh.np_elm,
                                self.vec_val, self.vec_velo)
    self.ls.Solve()
    self.vec_val += (self.ls.vec_x)*(self.dt*self.gamma_newmark) + (self.vec_velo)*self.dt
    self.vec_velo += self.ls.vec_x

  def step_time(self):
    self.solve()


class FEM_NavierStorks2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    ndimval = 3
    self.ls = FEM_LinSys(np,ndimval,mesh.psup())
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.dt = 0.1
    self.gamma_newmark = 0.6

  def solve(self):
    self.ls.SetZero()
    mergeLinSys_navierStorks2D(self.ls.mat, self.ls.vec_f,
                                1.0, 1000.0, 0.0, 0.0,
                                self.dt, self.gamma_newmark,
                                self.mesh.np_pos, self.mesh.np_elm,
                                self.vec_val, self.vec_velo)
    self.ls.Solve(is_asymmetric=True)
    self.vec_val += (self.ls.vec_x)*(self.dt*self.gamma_newmark) + (self.vec_velo)*self.dt
    self.vec_velo += self.ls.vec_x

  def step_time(self):
    self.solve()


