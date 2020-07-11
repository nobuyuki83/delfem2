####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import numpy, os, math

from .c_core import CppMatrixSparse, PreconditionerILU
from .c_core import cppAddMasterSlavePattern, matrixSquareSparse_setPattern, \
  cppMatSparse_ScaleBlk_LeftRight, \
  cppMatSparse_ScaleBlkLen_LeftRight
from .c_core import cppPrecILU_SetPattern_ILUk
from .c_core import \
  linearSystem_setMasterSlave, \
  linsys_solve_pcg, \
  masterSlave_distributeValue, \
  linsys_solve_bicgstab
from .c_core import \
  cppFEM_Merge_PointMass, \
  cppFEM_Merge_PointContact, \
  cppFEM_Merge_ScalarPoission, \
  cppFEM_Merge_ScalarDiffuse, \
  cppFEM_Merge_SolidLinearStatic, \
  cppFEM_Merge_SolidLinearDynamic, \
  cppFEM_Merge_FluidStorksStatic, \
  cppFEM_Merge_FluidStorksDynamic, \
  cppFEM_Merge_FluidNavierStorks, \
  cppFEM_Merge_ShellCloth, \
  cppFEM_Merge_ShellMitc3Static
from .c_core import \
  pbd_proj_rigid2d, \
  pbd_proj_rigid3d, \
  pbd_pointFixBC, \
  pbd_proj_cloth_stretch, \
  pbd_proj_cloth_bend, \
  pbd_proj_seam, \
  pbd_proj_contact
from .c_core import \
  write_vtk_meshpoint, \
  write_vtk_meshelem, \
  write_vtk_pointscalar, \
  write_vtk_pointvector  
from .c_core import matrixSquareSparse_setFixBC
from .c_core import elemQuad_dihedralTri, cppJArray_MeshPsup, jarray_add_diagonal, cppJArray_Sort
from .c_core import map_value
from .c_core import MathExpressionEvaluator
from .c_core import cppMassPoint_Mesh, cppMassLumped_ShellPlateBendingMitc3
from .c_core import CppSDF3

from .cadmsh import Mesh, MeshDynTri2D


def normalize_rigmsh(rigmsh):
  aabb = rigmsh.aabb()
  center = aabb.center()
  rigmsh.translate(-center[0],-center[1],-center[2])
  rigmsh.scale(1.0/aabb.max_length())

#####################################################

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
  def __init__(self):
    self.nlev_fill = 0
    self.nitr = 1000
    self.conv_ratio = 1.0e-4
    self.conv_hist = []

  def set_dimension(self,
                    np:int, ndimval:int):
    self.np = np
    self.ndimval = ndimval
    self.f = numpy.zeros((np,ndimval), dtype=numpy.float64)
    self.x = numpy.zeros((np,ndimval), dtype=numpy.float64)
    self.bc = numpy.zeros((np,ndimval), dtype=numpy.int32)
    self.ms = None
    self.mat = None
    self.mat_prec = None

  def set_pattern(self, pattern:tuple, master_slave_ptn=None):
    self.vec_ms = master_slave_ptn
    # matrix
    self.mat = CppMatrixSparse()
    self.mat.initialize(self.np, self.ndimval, True)
    psup_ind,psup = pattern[0],pattern[1]
    if self.vec_ms is not None:
      psup_ind1,psup1 = cppAddMasterSlavePattern(self.vec_ms,psup_ind,psup)
    else:
      psup_ind1,psup1 = psup_ind,psup
    assert psup_ind1.dtype == numpy.uint32 and psup1.dtype == numpy.uint32
    cppJArray_Sort(psup_ind1, psup1)
    matrixSquareSparse_setPattern(self.mat, psup_ind1, psup1)
    # preconditioner
    self.mat_prec = PreconditionerILU()
    cppPrecILU_SetPattern_ILUk(self.mat_prec, self.mat, self.nlev_fill)
#    cppPrecILU_SetPattern_ILUk(self.mat_prec, self.mat, 0)

  def set_zero(self):
    self.mat.set_zero()
    self.f[:,:] = 0.0

  def set_bc_ms(self):
    #### setting bc
    self.f[self.bc != 0] = 0.0
    matrixSquareSparse_setFixBC(self.mat, self.bc)
    if self.vec_ms is not None:
      linearSystem_setMasterSlave(self.mat,self.f, self.vec_ms)

  def set_precond(self):
    self.mat_prec.set_value(self.mat)
    self.mat_prec.ilu_decomp()

  def solve_iteration(self,is_asymmetric=False):
    #### solving matrix
    if not is_asymmetric:
      self.conv_hist = linsys_solve_pcg(self.f, self.x,
                                        self.conv_ratio, self.nitr,
                                        self.mat, self.mat_prec)
    else:
      self.conv_hist = linsys_solve_bicgstab(self.f, self.x,
                                             self.conv_ratio, self.nitr,
                                             self.mat, self.mat_prec)
    self.x[self.bc != 0] = 0.0

##########################################################################



class FEM_ScalarPoisson():
  def __init__(self,
               source=0.0,
               alpha=1.0):
    self.ls = FEM_LinSys()
    self.param_source = source
    self.param_alpha = alpha

  def updated_topology(self,mesh:Mesh,mapper=None,master_slave_pattern=None):
    self.mesh = mesh
    np = self.mesh.np_pos.shape[0]
    ndimval = 1
    val_new = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    if mapper is not None and hasattr(self,"value"):
      val_new[:self.value.shape[0],:] = self.value
      map_value(val_new,mapper)
    self.value = val_new
    self.ls.set_dimension(np,ndimval)
    self.ls.set_pattern(self.mesh.psup(),master_slave_ptn=master_slave_pattern)

  def solve(self):
    assert self.ls.mat is not None
    self.ls.set_zero()
    cppFEM_Merge_ScalarPoission(self.ls.mat, self.ls.f,
                                 self.param_alpha, self.param_source,
                                 self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                                 self.value)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration()
    self.value += self.ls.x
    if self.ls.vec_ms is not None:
      masterSlave_distributeValue(self.value, self.ls.vec_ms)


class FEM_ScalarDiffuse():
  def __init__(self):
    self.ls = FEM_LinSys()
    self.dt = 0.01
    self.gamma_newmark = 0.6
    self.param_alpha = 1.0
    self.param_rho = 1.0
    self.param_source = 1.0

  def updated_topology(self,mesh:Mesh):
    self.mesh = mesh
    np = self.mesh.np_pos.shape[0]
    ndimval = 1
    self.value = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.velocity = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls.set_dimension(np,ndimval)
    self.ls.set_pattern(self.mesh.psup())

  def step_time(self):
    assert self.ls.mat is not None
    self.ls.set_zero()
    cppFEM_Merge_ScalarDiffuse(self.ls.mat, self.ls.f,
                                self.param_alpha, self.param_rho, self.param_source,
                                self.dt, self.gamma_newmark,
                                self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                                self.value, self.velocity)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration()
    self.value += (self.ls.x)*(self.dt*self.gamma_newmark) + (self.velocity)*self.dt
    self.velocity += self.ls.x

  def initialize(self):
    self.value[:] = 0.0
    self.velocity[:] = 0.0

############################################################
## FEM solid from here
class FEM_SolidLinearStatic():
  def __init__(self):
    self.ls = FEM_LinSys()
    self.param_gravity_x = +0.0
    self.param_gravity_y = -0.0
    self.param_gravity_z = +0.0
    self.param_myu = 1.0
    self.param_lambda = 0.0
    self.param_rho = 1.0

  def updated_topology(self,mesh:Mesh):
    self.mesh = mesh
    np = self.mesh.np_pos.shape[0]
    ndimval = self.mesh.np_pos.shape[1]
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls.set_dimension(np, ndimval)
    self.ls.set_pattern(self.mesh.psup())

  def solve(self):
    assert self.ls.mat is not None
    self.ls.set_zero()
    gravity = [self.param_gravity_x, self.param_gravity_y]
    if self.mesh.np_pos.shape[1] == 3:
      gravity = [self.param_gravity_x, self.param_gravity_y, self.param_gravity_z]
    cppFEM_Merge_SolidLinearStatic(self.ls.mat, self.ls.f,
                                    self.param_myu, self.param_lambda, self.param_rho, gravity,
                                    self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                                    self.vec_val)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration()
    self.vec_val += self.ls.x


class FEM_SolidLinearEigen():
  def __init__(self):
    self.ls = FEM_LinSys()
    self.rho = 1.0
    self.mesh = None

  def updated_topology(self,mesh:Mesh):
    self.mesh = mesh
    np = self.mesh.np_pos.shape[0]
    ndimval = self.mesh.np_pos.shape[1]
    self.mode = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ker = numpy.zeros((6,np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.mass_lumped_sqrt_inv = numpy.zeros((np,), dtype=numpy.float64)
    self.ls.set_dimension(np, ndimval)
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.updated_geometry()

  def updated_geometry(self):
    cppMassPoint_Mesh(self.mass_lumped_sqrt_inv,
                       self.rho,
                       self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type)
    self.mass_lumped_sqrt_inv = numpy.sqrt(self.mass_lumped_sqrt_inv)
    if self.mesh.np_pos.shape[1] == 3:
      self.ker = self.ker.reshape((6,-1,3))
      self.ker[:,:,:] = 0.0
      self.ker[0,:,0] = + self.mass_lumped_sqrt_inv[:]
      self.ker[1,:,1] = + self.mass_lumped_sqrt_inv[:]
      self.ker[2,:,2] = + self.mass_lumped_sqrt_inv[:]
      self.ker[3,:,2] = - self.mass_lumped_sqrt_inv*self.mesh.np_pos[:,1]
      self.ker[3,:,1] = + self.mass_lumped_sqrt_inv*self.mesh.np_pos[:,2]
      self.ker[4,:,0] = - self.mass_lumped_sqrt_inv*self.mesh.np_pos[:,2]
      self.ker[4,:,2] = + self.mass_lumped_sqrt_inv*self.mesh.np_pos[:,0]
      self.ker[5,:,1] = - self.mass_lumped_sqrt_inv*self.mesh.np_pos[:,0]
      self.ker[5,:,0] = + self.mass_lumped_sqrt_inv*self.mesh.np_pos[:,1]
      self.ker = self.ker.reshape((6,-1))
      for i in range(6):
        self.ker[i] /= numpy.linalg.norm(self.ker[i])
        for j in range(i+1,6):
          self.ker[j] -= numpy.dot(self.ker[i], self.ker[j]) * self.ker[i]
      '''
      self.ker[0] /= numpy.linalg.norm(self.ker[0])
      self.ker[1] -= numpy.dot(self.ker[0],self.ker[1]) * self.ker[0]
      self.ker[2] -= numpy.dot(self.ker[0],self.ker[2]) * self.ker[0]
      self.ker[3] -= numpy.dot(self.ker[0],self.ker[3]) * self.ker[0]
      self.ker[4] -= numpy.dot(self.ker[0],self.ker[4]) * self.ker[0]
      self.ker[5] -= numpy.dot(self.ker[0],self.ker[5]) * self.ker[0]
      self.ker[1] /= numpy.linalg.norm(self.ker[1])
      self.ker[2] -= numpy.dot(self.ker[1],self.ker[2]) * self.ker[1]
      self.ker[3] -= numpy.dot(self.ker[1],self.ker[3]) * self.ker[1]
      self.ker[4] -= numpy.dot(self.ker[1],self.ker[4]) * self.ker[1]
      self.ker[5] -= numpy.dot(self.ker[1],self.ker[5]) * self.ker[1]
      self.ker[2] /= numpy.linalg.norm(self.ker[2])
      self.ker[3] -= numpy.dot(self.ker[2],self.ker[3]) * self.ker[2]
      self.ker[4] -= numpy.dot(self.ker[2],self.ker[4]) * self.ker[2]
      self.ker[5] -= numpy.dot(self.ker[2],self.ker[5]) * self.ker[2]
      self.ker[3] /= numpy.linalg.norm(self.ker[3])
      self.ker[4] -= numpy.dot(self.ker[3],self.ker[4]) * self.ker[3]
      self.ker[5] -= numpy.dot(self.ker[3],self.ker[5]) * self.ker[3]
      self.ker[4] /= numpy.linalg.norm(self.ker[4])
      self.ker[5] -= numpy.dot(self.ker[4],self.ker[5]) * self.ker[4]
      self.ker[5] /= numpy.linalg.norm(self.ker[5])
      '''
      '''
      for i in range(6):
        for j in range(i,6):
          print(i,j,numpy.dot(self.ker[i],self.ker[j]))
      '''
      self.ker = self.ker.reshape((6,-1,3))
      self.mass_lumped_sqrt_inv = numpy.reciprocal( self.mass_lumped_sqrt_inv )
    else:
      assert 0
      pass
    self.ls.set_zero()
    self.mode[:] = 0.0
    cppFEM_Merge_SolidLinearStatic(self.ls.mat, self.ls.f,
                                    1.0, 0.1, 0.0, (0,0,0),
                                    self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                                    self.mode)
    cppMatSparse_ScaleBlk_LeftRight(self.ls.mat, self.mass_lumped_sqrt_inv)
    self.ls.mat.add_dia(0.0)
    ####
    self.ls.set_precond()

  def solve(self):
    self.mode[:] = self.ls.f
    self.ls.x[:] = 0.0
    self.ls.solve_iteration()
    ####
    x = self.ls.x.reshape((-1))
    self.ker = self.ker.reshape((6,-1))
    x -= numpy.dot(x,self.ker[0]) * self.ker[0]
    x -= numpy.dot(x,self.ker[1]) * self.ker[1]
    x -= numpy.dot(x,self.ker[2]) * self.ker[2]
    x -= numpy.dot(x,self.ker[3]) * self.ker[3]
    x -= numpy.dot(x,self.ker[4]) * self.ker[4]
    x -= numpy.dot(x,self.ker[5]) * self.ker[5]
    x /= numpy.linalg.norm(x)
    self.ker = self.ker.reshape((6,-1,3))
    self.ls.f[:] = self.ls.x
    ####
    self.mode[:,0] = self.mass_lumped_sqrt_inv*self.ls.x[:,0]*0.03
    self.mode[:,1] = self.mass_lumped_sqrt_inv*self.ls.x[:,1]*0.03
    self.mode[:,2] = self.mass_lumped_sqrt_inv*self.ls.x[:,2]*0.03

  def step_time(self):
    self.solve()

  def initialize(self):
    self.ls.f[:] = numpy.random.uniform(-1, 1, self.mesh.np_pos.shape)


class FEM_SolidLinearDynamic():
  def __init__(self):
    self.ls = FEM_LinSys()
    self.mesh = None
    self.param_gravity_x = +0.0
    self.param_gravity_y = -0.0
    self.param_gravity_z = +0.0
    self.dt = 0.1
    self.gamma_newmark = 0.6
    self.beta_newmark = 0.36

  def updated_topology(self,mesh:Mesh):
    self.mesh = mesh
    np = self.mesh.np_pos.shape[0]
    ndimval = self.mesh.np_pos.shape[1]
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_acc = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls.set_dimension(np, ndimval)
    self.ls.set_pattern(self.mesh.psup())

  def solve(self):
    assert self.ls.mat is not None
    self.ls.set_zero()
    gravity = [self.param_gravity_x, self.param_gravity_y]
    if self.mesh.np_pos.shape[1] == 3:
      gravity = [self.param_gravity_x, self.param_gravity_y, self.param_gravity_z]
    cppFEM_Merge_SolidLinearDynamic(self.ls.mat, self.ls.f,
                                     1.0, 0.0, 1.0, gravity,
                                     self.dt, self.gamma_newmark, self.beta_newmark,
                                     self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                                     self.vec_val, self.vec_velo, self.vec_acc)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration()
    self.vec_val += (self.dt)*self.vec_velo + (0.5*self.dt*self.dt)*self.vec_acc + (self.dt*self.dt*self.beta_newmark)*self.ls.x
    self.vec_velo += (self.dt*self.gamma_newmark)*self.ls.x + (self.dt)*self.vec_acc
    self.vec_acc += self.ls.x

  def step_time(self):
    self.solve()


#############################################
## shell from here

class FEM_ShellPlateBendingMITC3():
  def __init__(self):
    self.ls = FEM_LinSys()
    self.mesh = None
    self.param_gravity_z = 0.001
    self.param_thickness = 0.05
    self.param_rho = 1.0
    self.param_myu = 100.0
    self.param_lambda = 100.0

  def updated_topology(self, mesh: Mesh):
    assert mesh.np_pos.shape[1] == 2
    self.mesh = mesh
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.disp = numpy.zeros((np,ndimval),dtype=numpy.float64)
    self.ls.set_dimension(np, ndimval)
    self.ls.set_pattern(self.mesh.psup())

  def solve(self):
    self.ls.set_zero()
    cppFEM_Merge_ShellMitc3Static(self.ls.mat, self.ls.f,
                                   self.param_thickness, self.param_lambda, self.param_myu,
                                   self.param_rho, self.param_gravity_z,
                                   self.mesh.np_pos, self.mesh.np_elm,
                                   self.disp)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration()
    self.disp += self.ls.x

  def initialize(self):
    self.disp[:] = 0.0


class FEM_ShellPlateBendingMITC3_Eigen():
  def __init__(self):
    self.ls = FEM_LinSys()
    self.param_thickness = 0.05
    self.param_rho = 1.0
    self.param_myu = 100.0
    self.param_lambda = 100.0
    self.param_offsetdia = 0.0
    self.freq_eigen = 0.0
    self.mesh = None

  def updated_topology(self, mesh: Mesh):
    self.mesh = mesh
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.mode = numpy.zeros((np, ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ker = numpy.zeros((3, np, ndimval), dtype=numpy.float64)  # initial guess is zero
    self.mass_lumped_sqrt_inv = numpy.zeros((np,3), dtype=numpy.float64)
    self.ls.set_dimension(np, ndimval)
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.updated_geometry()
    self.ls.f[:] = numpy.random.uniform(-1, 1, self.ls.f.shape)

  def updated_geometry(self):
    ####
    assert self.mass_lumped_sqrt_inv.shape[0] == self.mesh.np_pos.shape[0]
    cppMassLumped_ShellPlateBendingMitc3(self.mass_lumped_sqrt_inv,
                                         self.param_rho, self.param_thickness,
                                         self.mesh.np_pos,
                                         self.mesh.np_elm)
    self.mass_lumped_sqrt_inv = numpy.sqrt(self.mass_lumped_sqrt_inv)
    assert self.mesh.np_pos.shape[1] == 2
    self.ker = self.ker.reshape((3, -1, 3))
    self.ker[:, :, :] = 0.0
    self.ker[0, :, 0] = + self.mass_lumped_sqrt_inv[:,0]
    self.ker[1, :, 0] = + self.mass_lumped_sqrt_inv[:,0] * self.mesh.np_pos[:, 1]
    self.ker[1, :, 1] = + self.mass_lumped_sqrt_inv[:,1]
    self.ker[2, :, 0] = - self.mass_lumped_sqrt_inv[:,0] * self.mesh.np_pos[:, 0]
    self.ker[2, :, 2] = + self.mass_lumped_sqrt_inv[:,2]
    self.ker = self.ker.reshape((3, -1))
    for i in range(3):
      self.ker[i] /= numpy.linalg.norm(self.ker[i])
      for j in range(i + 1, 3):
        self.ker[j] -= numpy.dot(self.ker[i], self.ker[j]) * self.ker[i]
    self.ker = self.ker.reshape((3, -1, 3))
    self.mass_lumped_sqrt_inv = numpy.reciprocal(self.mass_lumped_sqrt_inv)
    self.ls.mat.set_zero()
    tmp = numpy.zeros_like(self.ls.f)
    cppFEM_Merge_ShellMitc3Static(self.ls.mat, tmp,
                                  self.param_thickness,
                                  self.param_lambda, self.param_myu,
                                  0.0, 0.0,
                                  self.mesh.np_pos, self.mesh.np_elm,
                                  self.mode)
    cppMatSparse_ScaleBlkLen_LeftRight(self.ls.mat,
                                       self.mass_lumped_sqrt_inv)
    self.ls.mat.add_dia(self.param_offsetdia)
    ###
    self.ls.set_precond()

  def solve(self):
    self.mode[:] = self.ls.f
    self.ls.x[:] = 0.0
    self.ls.solve_iteration()
    ####
    lam0 = numpy.dot(self.ls.x.flatten(),self.mode.flatten())
    if lam0 > 0.0 and 1.0 / lam0 - self.param_offsetdia > 0 :
      self.freq_eigen = math.sqrt(1.0/lam0-self.param_offsetdia)/(2.0*math.pi)
    else:
      self.freq_eigen = 0.0
    ####
    x = self.ls.x.reshape((-1))
    self.ker = self.ker.reshape((3, -1))
    x -= numpy.dot(x, self.ker[0]) * self.ker[0]
    x -= numpy.dot(x, self.ker[1]) * self.ker[1]
    x -= numpy.dot(x, self.ker[2]) * self.ker[2]
    nrm = numpy.linalg.norm(x)
    if abs(nrm) > 1.0e-10:
      x /= nrm
    self.ker = self.ker.reshape((3, -1, 3))
    self.ls.f[:] = self.ls.x
    ####
    self.mode[:] = self.mass_lumped_sqrt_inv[:] * self.ls.x[:] * 0.03

  def step_time(self):
    self.solve()

  def initialize(self):
    self.ls.f[:] = numpy.random.uniform(-1, 1, self.ls.f.shape)


class FEM_ShellCloth():
  def __init__(self):
    self.ls = FEM_LinSys()
    self.dt = 0.1
    self.gravity = (0,0,-1)
    self.rho = 1.0
    self.myu = 100.0
    self.lmd = 500.0
    self.sdf = None

  def updated_topology(self,mesh:Mesh,mapper=None):
    self.mesh = mesh
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    vec_val_new = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    vec_velo_new = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    if mapper is not None:
      vec_val_new[:self.vec_val.shape[0],:] = self.vec_val
      vec_velo_new[:self.vec_velo.shape[0],:] = self.vec_velo
      map_value(vec_val_new,mapper)
      map_value(vec_velo_new,mapper)
    else:
      vec_val_new[:,:2] = self.mesh.np_pos
    self.vec_val = vec_val_new
    self.vec_velo = vec_velo_new
    self.np_quad = elemQuad_dihedralTri(self.mesh.np_elm, np)
    self.ls.set_dimension(np,ndimval)
    self.ls.set_pattern(cppJArray_MeshPsup(self.np_quad, np))

  def solve(self):
    assert self.ls.mat is not None
    self.ls.set_zero()
    cppFEM_Merge_ShellCloth(self.ls.mat, self.ls.f,
                             self.myu, self.lmd, self.dt,
                             self.mesh.np_pos, self.mesh.np_elm,
                             self.np_quad,
                             self.vec_val)
    cppFEM_Merge_PointMass(self.ls.mat, self.ls.f,
                            self.rho, self.dt,
                            self.gravity,
                            self.vec_val, self.vec_velo)
    if self.sdf is not None:
      cppFEM_Merge_PointContact(self.ls.mat, self.ls.f,
                                 10000, 0.1,
                                 [self.sdf],
                                 self.vec_val)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration()
    self.vec_val += self.ls.x
    self.vec_velo = (1.0/self.dt)*self.ls.x

  def step_time(self):
    self.solve()


##########################
## fluid from here

class FEM_FluidStorksStatic():
  def __init__(self,
               mesh: Mesh):
    self.ls = FEM_LinSys()
    self.mesh = mesh
    self.updated_topology()

  def updated_topology(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls.set_dimension(np,ndimval)
    self.ls.set_pattern(self.mesh.psup())

  def solve(self):
    assert self.ls.mat is not None
    self.ls.set_zero()
    cppFEM_Merge_FluidStorksStatic(self.ls.mat, self.ls.f,
                                    1.0, 0.0, 0.0,
                                    self.mesh.np_pos, self.mesh.np_elm,
                                    self.vec_val)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration()
    self.vec_val += self.ls.x

  def step_time(self):
    self.solve()


class FEM_FluidStorksDynamic():
  def __init__(self,
               mesh: Mesh):
    self.ls = FEM_LinSys()
    self.dt = 0.005
    self.gamma_newmark = 0.6
    self.mesh = mesh
    self.updated_topology()

  def updated_topology(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls.set_dimension(np,ndimval)
    self.ls.set_pattern(self.mesh.psup())

  def solve(self):
    assert self.ls.mat is not None
    self.ls.set_zero()
    cppFEM_Merge_FluidStorksDynamic(self.ls.mat, self.ls.f,
                                     1.0, 1.0, 0.0, 0.0,
                                     self.dt, self.gamma_newmark,
                                     self.mesh.np_pos, self.mesh.np_elm,
                                     self.vec_val, self.vec_velo)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration()
    self.vec_val += (self.ls.x)*(self.dt*self.gamma_newmark) + (self.vec_velo)*self.dt
    self.vec_velo += self.ls.x

  def step_time(self):
    self.solve()


class FEM_FluidNavierStorks():
  def __init__(self,
               mesh: Mesh):
    self.ls = FEM_LinSys()
    self.dt = 0.1
    self.gamma_newmark = 0.6
    self.mesh = mesh
    self.updated_topology()

  def updated_topology(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls.set_dimension(np,ndimval)
    self.ls.set_pattern(self.mesh.psup())

  def solve(self):
    assert self.ls.mat is not None
    self.ls.set_zero()
    cppFEM_Merge_FluidNavierStorks(self.ls.mat, self.ls.f,
                                    1.0, 1000.0, 0.0, 0.0,
                                    self.dt, self.gamma_newmark,
                                    self.mesh.np_pos, self.mesh.np_elm,
                                    self.vec_val, self.vec_velo)
    self.ls.set_bc_ms()
    self.ls.set_precond()
    self.ls.solve_iteration(is_asymmetric=True)
    self.vec_val += (self.ls.x)*(self.dt*self.gamma_newmark) + (self.vec_velo)*self.dt
    self.vec_velo += self.ls.x

  def step_time(self):
    self.solve()

###########################################################################

class PBD():
  def __init__(self):
    self.dt = 0.1

  def updated_topology(self,mesh:Mesh):
    self.mesh = mesh
    np = mesh.np_pos.shape[0]
    self.vec_bc = numpy.zeros((np,), dtype=numpy.int32)
    self.vec_val = mesh.np_pos.copy()
    self.vec_velo = numpy.zeros_like(self.vec_val, dtype=numpy.float64)
    self.vec_tpos = mesh.np_pos.copy()
    self.psup = mesh.psup()
    self.psup = jarray_add_diagonal(*self.psup)

  def step_time(self):
    self.vec_tpos[:] = self.vec_val + self.dt * self.vec_velo
    pbd_pointFixBC(self.vec_tpos, self.vec_bc, self.vec_val)
    for itr in range(1):
      if self.mesh.np_pos.shape[1] == 2:
        pbd_proj_rigid2d(self.vec_tpos,
                     0.5, self.psup[0], self.psup[1],
                     self.mesh.np_pos)
      if self.mesh.np_pos.shape[1] == 3:
        pbd_proj_rigid3d(self.vec_tpos,
                     0.5, self.psup[0], self.psup[1],
                     self.mesh.np_pos)
    pbd_pointFixBC(self.vec_tpos, self.vec_bc, self.vec_val)
    self.vec_velo[:] = (self.vec_tpos-self.vec_val)/self.dt
    self.vec_val[:] = self.vec_tpos

  def initialize(self):
    self.vec_val[:] = self.mesh.np_pos[:]
    self.vec_velo[:] = 0.0


class PBD_Cloth():
  def __init__(self):
    self.dt = 0.1
    self.param_gravity_x = 0.0
    self.param_gravity_y = 0.0
    self.param_gravity_z = 0.0
    self.dmsh = None
    self.elems_seam = None
    self.sdf = None
    self.bc = None
    self.vec_val = None
    self.vec_velo = None
    self.vec_tpos = None

  def updated_topology(self, dmsh:MeshDynTri2D):
    self.dmsh = dmsh
    np = dmsh.np_pos.shape[0]
    self.bc = numpy.zeros((np,), dtype=numpy.int32)
    self.vec_val = numpy.zeros((np,3), dtype=numpy.float64)
    self.vec_val[:,:2] = dmsh.np_pos
    self.vec_velo = numpy.zeros_like(self.vec_val, dtype=numpy.float64)
    self.vec_tpos = self.vec_val.copy()

  def step_time(self):
    self.vec_tpos[:] = self.vec_val + self.dt * self.vec_velo
    self.vec_tpos[:,0] += self.dt*self.dt*self.param_gravity_x
    self.vec_tpos[:,1] += self.dt*self.dt*self.param_gravity_y
    self.vec_tpos[:,2] += self.dt*self.dt*self.param_gravity_z
    # pre
    pbd_pointFixBC(self.vec_tpos, self.bc, self.vec_val)
    # cloth
    for itr in range(1):
      pbd_proj_cloth_stretch(self.vec_tpos,self.dmsh.cdmsh)
    # bend
    for itr in range(1):
      pbd_proj_cloth_bend(self.vec_tpos, self.dmsh.cdmsh)
    # seam
    if isinstance(self.elems_seam,numpy.ndarray):
      for itr in range(1):
        pbd_proj_seam(self.vec_tpos,self.elems_seam)
    # contact
    if self.sdf is not None:
      self.sdf.project_points_outside(self.vec_tpos)
    # post
    pbd_pointFixBC(self.vec_tpos, self.bc, self.vec_val)
    self.vec_velo[:] = (self.vec_tpos-self.vec_val)/self.dt
    self.vec_val[:] = self.vec_tpos

  def initialize(self):
    self.vec_val[:,:2] = self.dmsh.np_pos
    self.vec_val[:,2] = 0.0
    self.vec_velo[:] = 0.0