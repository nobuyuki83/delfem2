####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import numpy, os
from typing import Tuple, List

from .c_core import CppCad2D, CppMeshDynTri2D, CppMesher_Cad2D, CppVoxelGrid, CppMapper, AABB3
from .c_core import CAD_EDGE_GEOM_LINE, CAD_EDGE_GEOM_BEZIER_CUBIC, CAD_EDGE_GEOM_BEZIER_QUADRATIC
from .c_core import cppCad2D_ImportSVG, cppSVG_Polyline
from .c_core import TRI, QUAD, HEX, TET, LINE
from .c_core import \
  meshquad2d_grid, \
  meshquad3d_voxelgrid, \
  meshquad3d_subdiv, \
  meshhex3d_voxelgrid, \
  meshhex3d_subdiv,\
  meshdyntri2d_initialize
from .c_core import meshtri3d_read_ply, meshtri3d_read_obj, meshtri3d_read_nastran, meshtri3d_write_obj
from .c_core import mvc
from .c_core import setXY_MeshDynTri2D
from .c_core import cad_getPointsEdge, cppJArray_MeshPsup, quality_meshTri2D
from .c_core import copyMeshDynTri2D
from .c_core import numpyXYTri_MeshDynTri2D
from .c_core import setTopology_ExtrudeTri2Tet
from .c_core import cppNormalVtx_Mesh, cppEdge_Mesh
from .c_core import \
  cppMeshTri3D_Cylinder, \
  cppMeshTri3D_Sphere, \
  cppMeshTri3D_GeoPoly, \
  cppMeshTri3D_Icosahedron, \
  cppMeshTri3D_Cube, \
  cppMeshTri3D_Torus
from .c_core import rotmat3_cartesian
from .c_core import num_node_elem

from .c_core import CppSDF3, CppSDF3_Sphere
from .c_core import project_points_outside_sdf
from .c_core import CppClliderPointsMeshTri3D

####################

class Mesh():
  def __init__(self,
               np_pos=numpy.ndarray((0,3),dtype=numpy.float64),
               np_elm=numpy.ndarray((0,3),dtype=numpy.uint32),
               elem_type=TRI):
    assert type(np_pos) == numpy.ndarray
    assert type(np_elm) == numpy.ndarray
    assert np_pos.dtype == numpy.float64
    assert np_elm.dtype == numpy.uint32
    assert np_elm.shape[1] == num_node_elem(elem_type)
    self.np_pos = np_pos
    self.np_elm = np_elm
    self.elem_type = elem_type

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

  #####

  def aabb3(self):
    return AABB3( self.minmax_xyz() )

  def ndim(self) -> int:
    return self.np_pos.shape[1]

  def scale_xyz(self,scale):
    self.np_pos *= scale

  def translate(self,d:List[float],mag=1.0):
    self.np_pos[:,0] += mag*d[0]
    self.np_pos[:,1] += mag*d[1]
    self.np_pos[:,2] += mag*d[2]

  def normalize(self, size=-1.0):
    bb = self.minmax_xyz()
    lenx = bb[1]-bb[0]
    leny = bb[3]-bb[2]
    lenz = bb[5]-bb[4]
    self.translate([-0.5*(bb[1]+bb[0]), -0.5*(bb[3]+bb[2]), -0.5*(bb[5]+bb[4])])
    lenmax = max(lenx,leny,lenz)
    if size < 0.0:
      size = 1.0
    self.scale_xyz(size/lenmax)

  def rotate(self,d:List[float]):
    R = rotmat3_cartesian(d)
    tmp0 = numpy.matmul(R, numpy.transpose(self.np_pos))
    self.np_pos[:,:] = numpy.transpose(tmp0)

  def subdiv(self):
    if self.elem_type == QUAD:
      np_pos1,np_quad1 = meshquad3d_subdiv(self.np_pos,self.np_elm)
      return Mesh(np_pos1,np_quad1,QUAD)
    if self.elem_type == HEX:
      np_pos1,np_quad1 = meshhex3d_subdiv(self.np_pos,self.np_elm)
      return Mesh(np_pos1,np_quad1,HEX)

  def psup(self) -> tuple:
    res = cppJArray_MeshPsup(self.np_elm, self.np_pos.shape[0])
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
    assert self.np_elm.dtype == numpy.uint32
    assert self.np_pos.dtype == numpy.float64

  def set_grid(self,shape:List[int]):
    if len(shape) == 2:
      self.np_pos, self.np_elm = meshquad2d_grid(shape[0], shape[1])
      self.elem_type = QUAD

  def set_extrude(self, msh0:'Mesh', nlayer:int):
    if msh0.elem_type ==TRI:
      self.elem_type = TET
      ####
      np0 = msh0.np_pos.shape[0]
      if self.np_pos.shape != (np0*(nlayer+1),3):
        self.np_pos = numpy.ndarray((np0*(nlayer+1),3),dtype=numpy.float64)
      for idiv in range(nlayer+1):
        self.np_pos[np0*idiv:np0*(idiv+1),:2] = msh0.np_pos[:,:]
        self.np_pos[np0*idiv:np0*(idiv+1), 2] = idiv
      ####
      ntri0 = msh0.np_elm.shape[0]
      if self.np_elm.shape != (ntri0*nlayer*3,4):
        self.np_elm = numpy.ndarray((ntri0*nlayer*3,4),dtype=numpy.uint32)
      setTopology_ExtrudeTri2Tet(self.np_elm,
                                 nlayer,np0,msh0.np_elm)

  def normal_vtx(self) -> numpy.ndarray:
    assert self.elem_type == TRI
    np_nrm = numpy.ndarray(self.np_pos.shape,dtype=numpy.float64)
    cppNormalVtx_Mesh(np_nrm,
                      self.np_pos, self.np_elm, self.elem_type)
    return np_nrm

  def mesh_edge(self):
    np_elem_edge = cppEdge_Mesh(self.np_pos, self.np_elm, self.elem_type)
    return Mesh(self.np_pos,np_elem_edge,LINE)

  def set_cylinder(self, r:float, l:float, nr:int, nl:int):
    self.np_pos, self.np_elm = cppMeshTri3D_Cylinder(r,l, nr, nl)
    self.elem_type = TRI

  def set_cube(self, n:int):
    self.np_pos, self.np_elm = cppMeshTri3D_Cube(n)
    self.elem_type = TRI

  def set_sphere(self, r:float, nla:int, nlo:int):
    self.np_pos, self.np_elm = cppMeshTri3D_Sphere(r, nla, nlo)
    self.elem_type = TRI

  def set_geopoly(self):
    self.np_pos, self.np_elm = cppMeshTri3D_GeoPoly()
    self.elem_type = TRI

  def set_icosahedron(self):
    self.np_pos, self.np_elm = cppMeshTri3D_Icosahedron()
    self.elem_type = TRI

  def set_torus(self, r0, r1):
    self.np_pos, self.np_elm = cppMeshTri3D_Torus(r0,r1)
    self.elem_type = TRI

###########################################################################

class SDF():
  def __init__(self):
    self.list_sdf = []

  def add(self,sdf):
    self.list_sdf.append(sdf)

  def project_points_outside(self,np_xyz):
    project_points_outside_sdf(np_xyz,self.list_sdf)

  def draw(self) -> None:
    for sdf in self.list_sdf:
      sdf.draw()


class Collider_PointsToMeshTri3D():
  def __init__(self):
    self.msh = None
    self.cpp_collider = CppClliderPointsMeshTri3D()

  def set_mesh(self,msh:Mesh):
    self.msh = msh
    self.cpp_collider.set_mesh(msh.np_pos,msh.np_elm,0.001)
    self.np_nrm = numpy.zeros_like(msh.np_pos,dtype=numpy.float64)
    cppNormalVtx_Mesh(self.np_nrm,
                      self.msh.np_pos, self.msh.np_elm, self.msh.elem_type)

  def project_points_outside(self,np_xyz):
    self.cpp_collider.project(np_xyz,
                              self.msh.np_pos, self.msh.np_elm,
                              self.np_nrm,
                              0.1)

###########################################################################


class MeshDynTri2D(Mesh):
  def __init__(self):
    super().__init__()
    self.elem_type = TRI
    self.np_pos = numpy.ndarray((0,2),dtype=numpy.float64)
    self.np_elm = numpy.ndarray((0,3),dtype=numpy.uint32)
    self.cdmsh = CppMeshDynTri2D()

  def set_mesh(self,msh:Mesh) -> None:
    assert msh.elem_type == TRI
    assert msh.np_pos.shape[1] == 2
    self.np_pos.resize(msh.np_pos.shape)
    self.np_elm.resize(msh.np_elm.shape)
    self.np_pos[:,:] = msh.np_pos
    self.np_elm[:,:] = msh.np_elm
    self.elem_type = msh.elem_type
    meshdyntri2d_initialize(self.cdmsh, self.np_pos, self.np_elm)

  def meshing_loops(self,loops:list,edge_length:float) -> None:
    self.cdmsh.meshing_loops(loops,edge_length)
    self.np_pos.resize((self.cdmsh.npoint(),2))
    self.np_elm.resize((self.cdmsh.ntri(),3))
    copyMeshDynTri2D(self.np_pos,self.np_elm, self.cdmsh)

  def refine_EdgeLongerThan_InsideCircle(self,elen,px,py,rad) -> CppMapper:
    mpr = CppMapper()
    self.cdmsh.refinementPlan_EdgeLongerThan_InsideCircle(mpr,
                                                          elen, px, py, rad)
    self.np_pos.resize((self.cdmsh.npoint(),2))
    self.np_elm.resize((self.cdmsh.ntri(),3))
    copyMeshDynTri2D(self.np_pos,self.np_elm, self.cdmsh)
    return mpr

  def syncXY_from_npPos(self) -> None:
    setXY_MeshDynTri2D(self.cdmsh, self.np_pos)

###########################################################################

class VoxelGrid:
  def __init__(self):
    self.vg = CppVoxelGrid()

  def add(self,ix,iy,iz):
    self.vg.add(ix,iy,iz, 1)

  def initialize(self, nx,ny,nz):
    self.vg.initialize(nx,ny,nz,0)

  def mesh_quad(self) -> Mesh:
    list_xyz, list_tri = meshquad3d_voxelgrid(self.vg)
    np_pos = numpy.array(list_xyz, dtype=numpy.float64).reshape((-1, 3))
    np_elm = numpy.array(list_tri, dtype=numpy.uint32).reshape((-1, 4))
    return Mesh(np_pos, np_elm, QUAD)

  def mesh_hex(self) -> Mesh:
    list_xyz, list_tri = meshhex3d_voxelgrid(self.vg)
    np_pos = numpy.array(list_xyz, dtype=numpy.float64).reshape((-1, 3))
    np_elm = numpy.array(list_tri, dtype=numpy.uint32).reshape((-1, 8))
    return Mesh(np_pos, np_elm, HEX)


###########################################################################


class Cad2D():
  def __init__(self):
    self.ccad = CppCad2D()

  def draw(self) -> None:
    self.ccad.draw()

  def mouse(self,btn,action,mods,src,dir,view_height) -> None:
    if btn == 0:
      if action == 1:
        self.ccad.pick(src[0],src[1],view_height)

  def motion(self,src0,src1,dir) -> None:
    self.ccad.drag_picked(src1[0],src1[1], src0[0],src0[1])

  def minmax_xyz(self):
    return self.ccad.minmax_xyz();

  ######

  def clear(self) -> None:
    """clear all the cad elements"""
    self.ccad.clear()

  def pick(self, x, y, view_height) -> None:
    self.ccad.pick(x,y,view_height)

  def add_polygon(self,list_xy) -> None:
    self.ccad.add_polygon(list_xy)
    self.ccad.check()

  def add_vtx_edge(self, iedge, pos:List[float]) -> None:
    self.ccad.add_vtx_edge(pos[0],pos[1],iedge)
    self.ccad.check()

  def add_vtx_face(self, iface, pos:List[float]) -> None:
    self.ccad.add_vtx_face(pos[0],pos[1],iface)
    self.ccad.check()

  def set_edge_type(self, iedge:int, type:int, param:List[float]):
    self.ccad.set_edge_type(iedge,type,param)

  def edge_type(self, iedge:int) -> int:
    return self.ccad.edge_type(iedge)

  def iedge_picked(self) -> int:
    return self.ccad.iedge_picked

  def ivtx_picked(self) -> int:
    return self.ccad.ivtx_picked

  def iface_picked(self) -> int:
    return self.ccad.iface_picked

  def clean_picked(self) -> None:
    self.ccad.ivtx_picked = -1
    self.ccad.iedge_picked = -1
    self.ccad.iface_picked = -1

  def points_edge(self, list_edge_index, np_xy, tolerance=0.01):
    return cad_getPointsEdge(self.ccad,list_edge_index, np_xy, tolerance=tolerance)

  def import_svg(self,path0:str,scale=(1.0,1.0)):
    if not os.path.isfile(path0):
      return false
    cppCad2D_ImportSVG(self.ccad,path0,scale[0],scale[1])
    self.ccad.check()

  def export_svg(self,path0:str,scale=1.0):
    list_xy = self.ccad.xy_vtxctrl_face(0)
    str0 = cppSVG_Polyline(list_xy,scale)
    with open(path0, mode='w') as f:
      f.write(str0)


########################################################################################

class CadMesh2D(Cad2D):
  def __init__(self,edge_length:float):
    super().__init__()
    self.ccad.is_draw_face = False
    self.edge_length = edge_length
    self.dmsh = MeshDynTri2D() # this object does not reallocate
    self.map_cad2msh = None # this object reallocate
    self.listW = list()
    self.is_sync_mesh = True
    self.mesher = Mesher_Cad2D(edge_length=edge_length)

  def draw(self):
    self.ccad.draw()
    self.dmsh.draw()

  def motion(self,src0,src1,dir):
    self.drag_picked(src1[0],src1[1], src0[0],src0[1])

  def minmax_xyz(self):
    return self.dmsh.minmax_xyz()

  #####

  def drag_picked(self, s1x,s1y, s0x,s0y):
    self.ccad.drag_picked(s1x,s1y, s0x,s0y)
    if not self.is_sync_mesh:
      return
    assert len(self.listW) == self.ccad.nface()
    for iface in range(self.ccad.nface()):
      list_xy_bound = self.ccad.xy_vtxctrl_face(iface)
      np_xy_bound = numpy.array(list_xy_bound).reshape([-1, 2])
      np_pos_face = numpy.dot(self.listW[iface][1],np_xy_bound)
      self.dmsh.np_pos[self.listW[iface][0]] = np_pos_face
      self.dmsh.syncXY_from_npPos()
      '''
      max_asp,min_area = quality_meshTri2D(self.dmsh.np_pos,self.dmsh.np_elm)
      if max_asp > 5.0 or min_area < 0.0:
        self.remesh()
      '''

  def remesh(self):
    self.mesher.meshing(self,self.dmsh)
    ####
    self.listW.clear()
    for iface in range(self.ccad.nface()):
      npIndPoint_face = self.mesher.points_on_faces([iface],self)
      npPosPoint_face = self.dmsh.np_pos[npIndPoint_face]
      np_xy_bound = numpy.array(self.ccad.xy_vtxctrl_face(iface)).reshape([-1, 2])
      W = mvc(npPosPoint_face, np_xy_bound)
      assert W.shape[0] == npPosPoint_face.shape[0]
      assert W.shape[1] == np_xy_bound.shape[0]
      self.listW.append( [npIndPoint_face,W] )
    assert len(self.listW) == self.ccad.nface()

  def add_vtx_edge(self, iedge:int, pos:List[float]):
    super().add_vtx_edge(iedge,[pos[0],pos[1]])
    self.remesh()

#  def add_polygon(self,list_xy):
#    self.ccad.add_polygon(list_xy)

#  def set_edge_type(self, iedge:int, type:int, param:List[float]):
#    super().set_edge_type(iedge,type,param)

#####################################################


class Mesher_Cad2D():
  def __init__(self,edge_length=0.01):
    self.cmshr = CppMesher_Cad2D()
    self.cmshr.edge_length = edge_length

  def points_on_faces(self,list_iface:List[int],cad:Cad2D) -> numpy.ndarray:
    list_points = self.cmshr.points_on_faces(list_iface,cad.ccad)
    return numpy.array(list_points,dtype=numpy.int32)

  def points_on_edges(self,list_iedge:List[int],cad:Cad2D) -> numpy.ndarray:
    list_points = self.cmshr.points_on_edges(list_iedge,cad.ccad)
    return numpy.array(list_points,dtype=numpy.int32)

  def points_on_one_edge(self,iedge:int,is_endpoints:bool, cad:Cad2D) -> numpy.ndarray:
    list_points = self.cmshr.points_on_one_edge(iedge,is_endpoints,cad.ccad)
    return numpy.array(list_points,dtype=numpy.int32)

  def meshing(self,cad:Cad2D,dmesh=None):
    cdmsh = CppMeshDynTri2D()
    self.cmshr.meshing(cdmsh,cad.ccad)
    np_pos, np_elm = numpyXYTri_MeshDynTri2D(cdmsh)
    if dmesh is None:
      dmesh = MeshDynTri2D()
    dmesh.cdmsh = cdmsh
    dmesh.np_pos = np_pos
    dmesh.np_elm = np_elm
    dmesh.elem_type = TRI
    return dmesh

