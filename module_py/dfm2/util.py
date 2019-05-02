import numpy, os
import OpenGL.GL as gl
from .dfm2 import *
from .cadmsh import *

def normalize_rigmsh(rigmsh):
  aabb = rigmsh.aabb()
  center = aabb.center()
  rigmsh.translate(-center[0],-center[1],-center[2])
  rigmsh.scale(1.0/aabb.max_length())

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
#      print(self.mesh.np_pos.shape,self.val_color.shape)
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


class FEM_Field():
  def __init__(self, fem,
               name_color="",
               name_disp=""):
    self.fem = fem
    ####
    self.name_color = name_color
    self.color_min = 0.0
    self.color_max = 0.3
    self.color_mode = 'bcgyr'
    ####
    self.name_disp = name_disp
    self.disp_mode = 'disp'

  def minmax_xyz(self):
    return self.fem.mesh.minmax_xyz()

  def draw(self):
    mesh = self.fem.mesh
    if hasattr(self.fem, self.name_color):
      npColor = getattr(self.fem, self.name_color)
      assert type(npColor) == numpy.ndarray
      self.color_map = ColorMap(self.color_min,self.color_max,self.color_mode)
      drawField_colorMap(mesh.np_pos, mesh.np_elm,
                         npColor,
                         self.color_map)

    if hasattr(self.fem, self.name_disp):
      npDisp = getattr(self.fem, self.name_disp)
      assert type(npDisp) == numpy.ndarray
      if self.disp_mode == 'disp':
        drawField_disp(mesh.np_pos, mesh.np_elm,
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


class FieldValueSetter():
  def __init__(self,
               mathexp: str,
               val:numpy.ndarray,
               idim:int,
               mesh:Mesh,
               npIdP: numpy.ndarray,
               dt:float):
    self.mesh = mesh
    self.val = val
    self.idim = idim
    self.mathexp = mathexp
    self.npIdP = npIdP
    self.eval = MathExpressionEvaluator()
    self.dt = dt
    self.time_cur = 0.0
    #####
    self.eval.set_key("x",0.0)
    self.eval.set_key("y",0.0)
    self.eval.set_key("z",0.0)
    self.eval.set_key("t",0.0)
    self.eval.set_expression(self.mathexp)

  def step_time(self):
    self.time_cur += self.dt
    for ip in self.npIdP:
      self.eval.set_key("t",self.time_cur)
      if self.mesh.np_pos.shape[1] >= 2:
        self.eval.set_key("x",self.mesh.np_pos[ip,0])
        self.eval.set_key("y",self.mesh.np_pos[ip,1])
      if self.mesh.np_pos.shape[1] == 3:
        self.eval.set_key("z",self.mesh.np_pos[ip,2])
      val = self.eval.eval()
      self.val[ip,self.idim] = self.mesh.np_pos[ip,self.idim] + val

######################################################

class FEM_LinSys():
  def __init__(self,
               np:int, ndimval:int):
    # vectors
    self.np = np
    self.ndimval = ndimval
    self.vec_f = numpy.zeros((np,ndimval), dtype=numpy.float64)
    self.vec_x = numpy.zeros((np,ndimval), dtype=numpy.float64)
    self.bc = numpy.zeros((np,ndimval), dtype=numpy.int32)
    self.vec_ms = numpy.ndarray((np,ndimval), dtype=numpy.int32)
    self.vec_ms.fill(-1)

    self.mat = None
    self.mat_prec = None

    self.nitr = 1000
    self.conv_ratio = 1.0e-4
    self.conv_hist = []

  def set_pattern(self, pattern:tuple):
    # matrix
    self.mat = MatrixSquareSparse()
    self.mat.initialize(self.np, self.ndimval, True)
    psup_ind,psup = pattern[0],pattern[1]
    psup_ind1,psup1 = addMasterSlavePattern(self.vec_ms,psup_ind,psup)
    jarray_sort(psup_ind1, psup1)
    matrixSquareSparse_setPattern(self.mat, psup_ind1, psup1)

    # preconditioner
    self.mat_prec = PreconditionerILU()
    precond_ilu0(self.mat_prec, self.mat)

  def set_zero(self):
    self.mat.setZero()
    self.vec_f[:,:] = 0.0

  def Solve(self,is_asymmetric=False):
    #### setting bc
    self.vec_f[self.bc != 0] = 0.0
    matrixSquareSparse_setFixBC(self.mat, self.bc)
    linearSystem_setMasterSlave(self.mat,self.vec_f, self.vec_ms)

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
    self.vec_x[self.bc != 0] = 0.0

##########################################################################



class FEM_Poisson():
  def __init__(self,
               mesh: Mesh,
               source=0.0):
    self.alpha = 1.0
    self.source = source
    self.mesh = mesh
    self.updated_mesh()

  def updated_mesh(self,mapper=None):
    np = self.mesh.np_pos.shape[0]
    ndimval = 1
    val_new = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    if mapper is not None:
      map_value(val_new,self.vec_val,mapper)
    self.value = val_new
    self.ls = FEM_LinSys(np,ndimval)

  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
    mergeLinSys_poission(self.ls.mat, self.ls.vec_f,
                         self.alpha, self.source,
                         self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                         self.value)
    self.ls.Solve()
    self.value += self.ls.vec_x
    masterSlave_distributeValue(self.value, self.ls.vec_ms)


class FEM_Diffuse():
  def __init__(self,
               mesh: Mesh,
               source=0.0):
    self.mesh = mesh
    self.dt = 0.01
    self.gamma_newmark = 0.6
    self.alpha = 1.0
    self.rho = 1.0
    self.source = source
    self.updated_mesh()

  def updated_mesh(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 1
    self.value = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.velocity = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls = FEM_LinSys(np,ndimval)

  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
    mergeLinSys_diffuse(self.ls.mat, self.ls.vec_f,
                        self.alpha, self.rho, self.source,
                        self.dt, self.gamma_newmark,
                        self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                        self.value, self.velocity)
    self.ls.Solve()
    self.value += (self.ls.vec_x)*(self.dt*self.gamma_newmark) + (self.velocity)*self.dt
    self.velocity += self.ls.vec_x

  def step_time(self):
    self.solve()


class FEM_LinearSolidStatic():
  def __init__(self,
               mesh: Mesh,
               gravity = [0,0,0]):
    self.mesh = mesh
    self.gravity = gravity
    self.updated_mesh()

  def updated_mesh(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = self.mesh.np_pos.shape[1]
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls = FEM_LinSys(np, ndimval)


  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
    mergeLinSys_linearSolidStatic(self.ls.mat, self.ls.vec_f,
                                  1.0, 0.0, 1.0, self.gravity,
                                  self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                                  self.vec_val)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x


class FEM_LinearSolidDynamic():
  def __init__(self,
               mesh: Mesh,
               gravity=[0, 0, 0]):
    self.mesh = mesh
    self.gravity = gravity
    self.dt = 0.1
    self.gamma_newmark = 0.6
    self.beta_newmark = 0.36
    self.updated_mesh()

  def updated_mesh(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = self.mesh.np_pos.shape[1]
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_acc = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls = FEM_LinSys(np, ndimval)

  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
    mergeLinSys_linearSolidDynamic(self.ls.mat, self.ls.vec_f,
                                   1.0, 0.0, 1.0, self.gravity,
                                   self.dt, self.gamma_newmark, self.beta_newmark,
                                   self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                                   self.vec_val, self.vec_velo, self.vec_acc)
    self.ls.Solve()
    self.vec_val += (self.dt)*self.vec_velo + (0.5*self.dt*self.dt)*self.vec_acc + (self.dt*self.dt*self.beta_newmark)*self.ls.vec_x
    self.vec_velo += (self.dt*self.gamma_newmark)*self.ls.vec_x + (self.dt)*self.vec_acc
    self.vec_acc += self.ls.vec_x

  def step_time(self):
    self.solve()


class FEM_Cloth():
  def __init__(self,mesh):
    self.dt = 0.1
    self.mesh = mesh
    self.sdf = SDF()
    self.updated_mesh()

  def updated_mesh(self,mapper=None):
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    vec_val_new = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    vec_velo_new = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    if mapper is not None:
      map_value(vec_val_new,self.vec_val,mapper)
      map_value(vec_velo_new,self.vec_velo,mapper)
    else:
      vec_val_new[:,:2] = self.mesh.np_pos
    self.vec_val = vec_val_new
    self.vec_velo = vec_velo_new
    self.ls = FEM_LinSys(np,ndimval)

  def solve(self):
    if self.ls.mat is None:
      np = self.mesh.np_pos.shape[0]
      self.np_quad = elemQuad_dihedralTri(self.mesh.np_elm, np)
      self.ls.set_pattern(jarray_mesh_psup(self.np_quad, np))
    self.ls.set_zero()
    mergeLinSys_cloth(self.ls.mat, self.ls.vec_f,
                      10.0, 500.0, self.dt,
                      self.mesh.np_pos, self.mesh.np_elm,
                      self.np_quad,
                      self.vec_val)
    mergeLinSys_massPoint(self.ls.mat, self.ls.vec_f,
                          1.0, self.dt,
                          [0,0,-1],
                          self.vec_val, self.vec_velo)
    mergeLinSys_contact(self.ls.mat, self.ls.vec_f,
                        10000, 0.1,
                        self.sdf.list_sdf,
                        self.vec_val)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x
    self.vec_velo = (1.0/self.dt)*self.ls.vec_x

  def step_time(self):
    self.solve()


##########################
## fluid from here

class FEM_StorksStatic2D():
  def __init__(self,
               mesh: Mesh):
    self.mesh = mesh
    self.updated_mesh()

  def updated_mesh(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.ls = FEM_LinSys(np,ndimval)
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero

  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
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
    self.dt = 0.005
    self.gamma_newmark = 0.6
    self.mesh = mesh
    self.updated_mesh()

  def updated_mesh(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.ls = FEM_LinSys(np,ndimval)
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero

  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
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
    self.dt = 0.1
    self.gamma_newmark = 0.6
    self.mesh = mesh
    self.updated_mesh()

  def updated_mesh(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls = FEM_LinSys(np,ndimval)

  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
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


class PBD():
  def __init__(self,
               mesh: Mesh):
    np = mesh.np_pos.shape[0]
    self.mesh = mesh
    self.vec_bc = numpy.zeros((np,), dtype=numpy.int32)
    self.vec_val = mesh.np_pos.copy()
    self.vec_velo = numpy.zeros_like(self.vec_val, dtype=numpy.float64)
    self.vec_tpos = mesh.np_pos.copy()
    self.dt = 0.1
    self.psup = mesh.psup()
    self.psup = jarray_add_diagonal(*self.psup)

  def step_time(self):
    self.vec_tpos[:] = self.vec_val + self.dt * self.vec_velo
    pointFixBC(self.vec_tpos, self.vec_bc, self.vec_val)
    for itr in range(1):
      if self.mesh.np_pos.shape[1] == 2:
        proj_rigid2d(self.vec_tpos,
                     0.5, self.psup[0], self.psup[1],
                     self.mesh.np_pos)
      if self.mesh.np_pos.shape[1] == 3:
        proj_rigid3d(self.vec_tpos,
                     0.5, self.psup[0], self.psup[1],
                     self.mesh.np_pos)
    pointFixBC(self.vec_tpos, self.vec_bc, self.vec_val)
    self.vec_velo[:] = (self.vec_tpos-self.vec_val)/self.dt
    self.vec_val[:] = self.vec_tpos

