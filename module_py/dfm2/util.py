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
               np_elm=numpy.ndarray((0,3),dtype=numpy.int32),
               elem_type=None):
    print("PyMesh -- construct")
    assert type(np_pos) == numpy.ndarray
    assert type(np_elm) == numpy.ndarray
    assert np_pos.dtype == numpy.float64
    assert np_elm.dtype == numpy.int32
    self.color_face = [0.8, 0.8, 0.8, 1.0]
    self.is_draw_edge = True
    self.np_pos = np_pos
    self.np_elm = np_elm
    self.elem_type = elem_type

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
    np_pos1,np_quad1 = meshquad3d_subdiv(self.np_pos,self.np_elm)
    return Mesh(np_pos1,np_quad1,QUAD)

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
      np_pos, np_elm = meshtri3d_read_ply(path_file)
    if ext == 'obj':
      np_pos, np_elm = meshtri3d_read_obj(path_file)
    if ext == 'nas' or ext == 'bdf':
      np_pos, np_elm = meshtri3d_read_nastran(path_file)
    return Mesh(np_pos,np_elm,TRI)
  return None

def mesh_voxelgrid(voxelgrid) -> Mesh:
  list_xyz, list_tri = getmesh_voxelgrid(voxelgrid)
  np_pos = numpy.array(list_xyz, dtype=numpy.float64).reshape((-1, 3))
  np_elm = numpy.array(list_tri, dtype=numpy.int32).reshape((-1, 4))
  return Mesh(np_pos,np_elm)

def mesh_grid(shape) -> Mesh:
  np_pos,np_quad = meshquad2d_grid(shape[0],shape[1])
  return Mesh(np_pos,np_quad,QUAD)

def mesh_cad(cad,len) -> Mesh:
  xy,tri = getMesh_cad(cad,len)
  mesh = Mesh(xy,tri,TRI)
  return mesh


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

  def write_vtk(self, path_vtk, message=""):
    write_vtk_meshpoint(path_vtk,"foobar", self.mesh.np_pos)
    write_vtk_meshelem(path_vtk, self.mesh.np_elm, self.mesh.elem_type)
    open(path_vtk, "a+").write("POINT_DATA {0}\n".format(self.val_color.shape[0]))
    print(self.val_color.shape)
    if self.val_color.ndim == 1:
      write_vtk_pointscalar(path_vtk, self.val_color)
    if self.val_color.ndim == 2 and self.val_color.shape[1] == 1:
      write_vtk_pointscalar(path_vtk, self.val_color)


class VisFEM_Color():
  def __init__(self):
    self.draw_val_min = 0.0
    self.draw_val_max = 1.0
    self.color_mode = 'bcgyr'
    self.is_update_min_max = True
    self.color_map = ColorMap(self.draw_val_min, self.draw_val_max, self.color_mode)

  def update(self,mesh,val_color):
    if self.is_update_min_max:
      self.draw_val_min = val_color.min()
      self.draw_val_max = val_color.max()
      self.color_map = ColorMap(self.draw_val_min, self.draw_val_max, self.color_mode)

  def draw(self,mesh,val_color):
    drawField_colorMap(mesh.np_pos, mesh.np_elm,
                       val_color,
                       self.color_map)

######################################################

class FEM_LinSys():
  def __init__(self,
               np:int, ndimval:int):
    # vectors
    self.np = np
    self.ndimval = ndimval
    self.vec_f = numpy.zeros((np,ndimval), dtype=numpy.float64)
    self.vec_x = numpy.zeros((np,ndimval), dtype=numpy.float64)
    self.vec_bc = numpy.zeros((np,ndimval), dtype=numpy.int32)
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
    sortIndexedArray(psup_ind1, psup1)
    matrixSquareSparse_setPattern(self.mat, psup_ind1, psup1)

    # preconditioner
    self.mat_prec = PreconditionerILU()
    precond_ilu0(self.mat_prec, self.mat)

  def set_zero(self):
    self.mat.setZero()
    self.vec_f[:,:] = 0.0

  def Solve(self,is_asymmetric=False):
    #### setting bc
    self.vec_f[self.vec_bc != 0] = 0.0
    matrixSquareSparse_setFixBC(self.mat, self.vec_bc)
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
    self.vec_x[self.vec_bc != 0] = 0.0

##########################################################################



class FEM_Poisson():
  def __init__(self,
               mesh: Mesh,
               source=0.0):
    self.alpha = 1.0
    self.source = source
    self.mesh = mesh
    self.updated_mesh()

  def updated_mesh(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 1
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls = FEM_LinSys(np,ndimval)

  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
    mergeLinSys_poission(self.ls.mat, self.ls.vec_f,
                         self.alpha, self.source,
                         self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                         self.vec_val)
    self.ls.Solve()
    self.vec_val += self.ls.vec_x
    masterSlave_distributeValue(self.vec_val, self.ls.vec_ms)



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
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.ls = FEM_LinSys(np,ndimval)

  def solve(self):
    if self.ls.mat is None:
      self.ls.set_pattern(self.mesh.psup())
    self.ls.set_zero()
    mergeLinSys_diffuse(self.ls.mat, self.ls.vec_f,
                        self.alpha, self.rho, self.source,
                        self.dt, self.gamma_newmark,
                        self.mesh.np_pos, self.mesh.np_elm, self.mesh.elem_type,
                        self.vec_val, self.vec_velo)
    self.ls.Solve()
    self.vec_val += (self.ls.vec_x)*(self.dt*self.gamma_newmark) + (self.vec_velo)*self.dt
    self.vec_velo += self.ls.vec_x

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
    self.updated_mesh()

  def updated_mesh(self):
    np = self.mesh.np_pos.shape[0]
    ndimval = 3
    self.vec_val = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_velo = numpy.zeros((np,ndimval), dtype=numpy.float64)  # initial guess is zero
    self.vec_val[:,:2] = self.mesh.np_pos
    self.ls = FEM_LinSys(np,ndimval)

  def solve(self):
    if self.ls.mat is None:
      np = self.mesh.np_pos.shape[0]
      self.np_quad = elemQuad_dihedralTri(self.mesh.np_elm, np)
      self.ls.set_pattern(psup_mesh(self.np_quad, np))
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
    self.ls.Solve()
    self.vec_val += self.ls.vec_x
    self.vec_velo = (1.0/self.dt)*self.ls.vec_x

  def step_time(self):
    self.solve()


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


