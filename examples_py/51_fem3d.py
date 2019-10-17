####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import numpy
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def make_mesh():
  sdf0 = dfm2.CppSDF3_Sphere(0.55,[-0.5,0,0],True)
  sdf1 = dfm2.CppSDF3_Sphere(0.55,[+0.5,0,0],True)
  np_xyz,np_tet = dfm2.isosurface([sdf0,sdf1])
  print(np_xyz.shape,np_tet.shape)
  msh = dfm2.Mesh(np_xyz,np_tet,dfm2.TET)
  return msh

def poission(msh,npIdP0,npIdP1):
  fem = dfm2.FEM_ScalarPoisson()
  fem.updated_topology(msh)
  fem.ls.bc[npIdP0] = 1
  fem.ls.bc[npIdP1] = 2
  fem.value[:] = 0.5
  fem.value[npIdP0] = 0.0
  fem.value[npIdP1] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  vis_color = dfm2.gl.VisFEM_ColorContour(fem,name_color="value")
  vis_color.set_color_minmax()
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([vis_color,axis])

def diffuse(msh,npIdP0,npIdP1):
  fem = dfm2.FEM_ScalarDiffuse()
  fem.updated_topology(msh)
  fem.ls.bc[npIdP1] = 1
  fem.value[:] = 0.0
  fem.value[npIdP1] = 1.0
  ####
  vis_color = dfm2.gl.VisFEM_ColorContour(fem,name_color="value")
  vis_color.draw_val_min = 0.0
  vis_color.draw_val_max = 1.0
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([fem,vis_color,axis])


def linear_solid_static(msh,npIdP):
  fem = dfm2.FEM_SolidLinearStatic()
  fem.param_gravity_x = +0.3
  fem.updated_topology(msh)
  fem.ls.bc[npIdP,:] = 1
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  vis_disp = dfm2.gl.VisFEM_ColorContour(fem,name_disp="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([vis_disp,axis])
  vis_disp.write_vtk("linearsolid3d.vtk")


def linear_solid_dynamic(msh,npIdP):
  fem = dfm2.FEM_SolidLinearDynamic()
  fem.param_gravity_x = +0.3
  fem.updated_topology(msh)
  fem.ls.bc[npIdP,:] = 1
  ####
  vis_disp = dfm2.gl.VisFEM_ColorContour(fem,name_disp="vec_val")
  axis = dfm2.gl.AxisXYZ(1.0)
  dfm2.gl.glfw.winDraw3d([fem,vis_disp,axis])


def main():
  msh = make_mesh()
  npIdP0 = numpy.where(msh.np_pos[:,0]>+1)
  npIdP1 = numpy.where(msh.np_pos[:,0]<-1)
  poission(msh,npIdP0,npIdP1)
  diffuse(msh,npIdP0,npIdP1)
  linear_solid_static(msh,npIdP1)
  linear_solid_dynamic(msh,npIdP1)


if __name__ == "__main__":
  main()