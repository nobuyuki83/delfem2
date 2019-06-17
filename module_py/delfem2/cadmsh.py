####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import numpy, os
import OpenGL.GL as gl
from typing import Tuple, List

from .libdelfem2 import CppCad2D
from .libdelfem2 import TRI, QUAD, HEX, TET
from .libdelfem2 import \
  meshquad2d_grid, \
  meshquad3d_voxelgrid, \
  meshquad3d_subdiv, \
  meshhex3d_voxelgrid, \
  meshhex3d_subdiv
from .libdelfem2 import draw_mesh_facenorm, draw_mesh_edge
from .libdelfem2 import meshtri3d_read_ply, mvc, meshtri3d_read_obj, meshdyntri2d_initialize
from .libdelfem2 import meshDynTri2D_CppCad2D, cad_getPointsEdge, jarray_mesh_psup, quality_meshTri2D
from .libdelfem2 import CppMeshDynTri2D, copyMeshDynTri2D, CppMapper
from .libdelfem2 import CppVoxelGrid, numpyXYTri_MeshDynTri2D

####################

class Mesh():
  def __init__(self,
               np_pos=numpy.ndarray((0,3),dtype=numpy.float64),
               np_elm=numpy.ndarray((0,3),dtype=numpy.int32),
               elem_type=TRI):
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

  def psup(self) -> tuple:
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

##########################################################################

def mesh_grid(shape:list) -> Mesh:
  np_pos,np_quad = meshquad2d_grid(shape[0],shape[1])
  return Mesh(np_pos,np_quad,QUAD)


###########################################################################


class MeshDynTri2D(Mesh):
  def __init__(self):
    super().__init__()
    self.elem_type = TRI
    self.np_pos = numpy.ndarray((0,2),dtype=numpy.float64)
    self.np_elm = numpy.ndarray((0,3),dtype=numpy.int32)
    self.cdmsh = CppMeshDynTri2D()

  def set_mesh(self,msh=Mesh) -> None:
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


###########################################################################


class Cad2D():
  def __init__(self):
    self.ccad = CppCad2D()

  def draw(self) -> None:
    self.ccad.draw()

  def mouse(self,btn,action,mods,src,dir,view_height):
    if btn == 0:
      if action == 1:
        self.ccad.pick(src[0],src[1],view_height)
      elif action == 0:
        self.clean_picked()

  def motion(self,src0,src1,dir) -> None:
    self.ccad.drag_picked(src1[0],src1[1], src0[0],src0[1])

  ########

  def pick(self, x, y, view_height) -> None:
    self.ccad.pick(x,y,view_height)

  def add_polygon(self,list_xy) -> None:
    self.ccad.add_polygon(list_xy)
    self.ccad.check()

  def add_vtx_edge(self, x, y, iedge) -> None:
    self.ccad.add_vtx_edge(x,y,iedge)
    self.ccad.check()

  def mesh(self,edge_len=0.05) -> Tuple[Mesh,List[numpy.ndarray]]:
    cdmsh, flg_pnt, flg_tri = meshDynTri2D_CppCad2D(self.ccad, edge_len)
    np_pos, np_elm = numpyXYTri_MeshDynTri2D(cdmsh)
    dmesh = MeshDynTri2D()
    dmesh.cdmsh = cdmsh
    dmesh.np_pos = np_pos
    dmesh.np_elm = np_elm
    dmesh.elem_type = TRI
    return dmesh, [flg_pnt, flg_tri]

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

######################

class CadMesh2D(Cad2D):
  def __init__(self,edge_length:float):
    super().__init__()
    self.ccad.is_draw_face = False
    self.edge_length = edge_length
    self.dmsh = MeshDynTri2D()
    self.map_cad2msh = None
    self.listW = list()

  def draw(self):
    self.ccad.draw()
    self.dmsh.draw()

  def motion(self,src0,src1,dir):
    self.drag_picked(src1[0],src1[1], src0[0],src0[1])

  def minmax_xyz(self):
    return self.msh.minmax_xyz()

  #####

  def drag_picked(self, s1x,s1y, s0x,s0y):
    self.ccad.drag_picked(s1x,s1y, s0x,s0y)
    assert len(self.listW) == self.ccad.nface()
    for iface in range(self.ccad.nface()):
      list_xy_bound = self.ccad.xy_vtx_face(iface)
      np_xy_bound = numpy.array(list_xy_bound).reshape([-1, 2])
      np_pos_face = numpy.dot(self.listW[iface][1],np_xy_bound)
      self.dmsh.np_pos[self.listW[iface][0]] = np_pos_face
    max_asp,min_area = quality_meshTri2D(self.dmsh.np_pos,self.dmsh.np_elm)
    if max_asp > 5.0 or min_area < 0.0:
      self.remesh()

  def remesh(self):
    self.dmsh.cdmsh, flg_pnt, flg_tri = meshDynTri2D_CppCad2D(self.ccad, self.edge_length)
    self.dmsh.np_pos, self.dmsh.np_elm = numpyXYTri_MeshDynTri2D(self.dmsh.cdmsh)
    self.map_cad2msh = [flg_pnt,flg_tri]
    ####
    self.listW.clear()
    for iface in range(self.ccad.nface()):
      np_ind_face = numpy.where(self.map_cad2msh[0] == self.ccad.nvtx()+self.ccad.nedge()+iface)[0]
      for ie,dir in self.ccad.ind_edge_face(iface):
        np_ind_face = numpy.append(np_ind_face, numpy.where( self.map_cad2msh[0] == self.ccad.nvtx()+ie )[0] )
      for iv in self.ccad.ind_vtx_face(iface):
        np_ind_face = numpy.append(np_ind_face, numpy.where( self.map_cad2msh[0] == iv )[0] )
      np_pos_face = self.dmsh.np_pos[np_ind_face]
      np_xy_bound = numpy.array(self.ccad.xy_vtx_face(iface)).reshape([-1, 2])
      W = mvc(np_pos_face, np_xy_bound)
      assert W.shape[0] == np_pos_face.shape[0]
      assert W.shape[1] == np_xy_bound.shape[0]
      self.listW.append( [np_ind_face,W] )
    assert len(self.listW) == self.ccad.nface()

  def add_vtx_edge(self, x, y, iedge):
    super().add_vtx_edge(x,y,iedge)
    self.remesh()

  def add_polygon(self,list_xy):
    self.ccad.add_polygon(list_xy)
    self.remesh()

#####################################################

class SDF():
  def __init__(self):
    self.list_sdf = []

  def draw(self):
    for sdf in self.list_sdf:
      sdf.draw()