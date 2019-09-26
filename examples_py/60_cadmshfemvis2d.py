####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw
import PyDelFEM2.cadmshsimvis
import numpy

def mesh():
  cmf = dfm2.CadMesh2D(edge_length=0.1)
  cmf.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  cmf.set_edge_type(0,1,[0.2, 0.3, -0.2, 0.3])
  cmf.remesh()
  dfm2.gl.glfw.winDraw3d([cmf])

def poisson():
  cmf = dfm2.cadmshsimvis.CadMesh2D_FEMPoisson(edge_length=0.1)
  cmf.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  cmf.remesh()
  dfm2.gl.glfw.winDraw3d([cmf])

def solid_linear_static():
  cmf = dfm2.cadmshsimvis.CadMesh2D_FEMSolidLinearStatic(edge_length=0.1)
  cmf.fem.param_gravity_y = -0.1
  cmf.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  cmf.remesh()
  dfm2.gl.glfw.winDraw3d([cmf])

def solid_linear_eigen():
  cmf = dfm2.cadmshsimvis.CadMesh2D_FEMSolidLinearEigen(edge_length=0.1)
  cmf.add_polygon([-1,-0.2, +1,-0.2, +1,+0.2, -1,+0.2])
  cmf.remesh()
  dfm2.gl.glfw.winDraw3d([cmf])

def shell_platebending_mitc3_eigen():
  cmf = dfm2.cadmshsimvis.CadMesh2D_FEMShellPlateBendingMITC3Eigen(edge_length=0.1)
  cmf.add_polygon([-1,-0.2, +1,-0.2, +1,+0.2, -1,+0.2])
  cmf.remesh()
  dfm2.gl.glfw.winDraw3d([cmf])


def pbd2d():
  cmf = dfm2.cadmshsimvis.CadMesh2D_PBD(edge_length=0.2)
  cmf.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  cmf.remesh()
  dfm2.gl.glfw.winDraw3d([cmf])

def pbd_cloth():
  cms = dfm2.cadmshsimvis.CadMesh2D_PBDCloth(edge_length=0.05)
  cms.add_polygon([0, 0, +1, 0, +1, +1, 0, +1])
  cms.add_polygon([2, 0, +3, 0, +3, +1, 2, +1])
  cms.list_seam = [[1,7],[3,5]]
  cms.remesh()

  trans0 = dfm2.Trans_Rigid2DTo3D(org2=[0.5, 0.5], org3=[0.0,0.0,0.5])
  npIndP_Face0 = cms.mesher.points_on_faces([0],cms)
  cms.pbd.vec_val[npIndP_Face0] = trans0.trans(cms.dmsh.np_pos[npIndP_Face0])

  trans1 = dfm2.Trans_Rigid2DTo3D(org2=[2.5,0.5], org3=[0.0,0.0,-0.5])
  trans1.R = dfm2.util.mat3_rot_cartesian(numpy.array([0,3.1415,0]))
  npIndP_Face1 = cms.mesher.points_on_faces([1],cms)
  cms.pbd.vec_val[npIndP_Face1] = trans1.trans(cms.dmsh.np_pos[npIndP_Face1])

  cms.pbd.sdf = dfm2.SDF()
  cms.pbd.sdf.add( dfm2.CppSDF3_Sphere(0.3,[0,0,0],True) )

  dfm2.gl.glfw.winDraw3d([cms,cms.pbd.sdf])


if __name__ == "__main__":
  shell_platebending_mitc3_eigen()
  mesh()
  poisson()
  solid_linear_static()
  solid_linear_eigen()
  pbd2d()
  pbd_cloth()