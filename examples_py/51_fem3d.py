####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import sys, numpy
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw

def make_mesh():
  sdf = dfm2.SDF()
  sdf.list_sdf.append( dfm2.SDF_Sphere(0.55,[-0.5,0,0],True) )
  sdf.list_sdf.append( dfm2.SDF_Sphere(0.55,[+0.5,0,0],True) )
  np_xyz,np_tet = dfm2.isosurface(sdf.list_sdf)
  print(np_xyz.shape,np_tet.shape)
  msh = dfm2.Mesh(np_xyz,np_tet,dfm2.TET)
  return msh

def poission(msh,npIdP0,npIdP1):
  fem = dfm2.FEM_Poisson(msh)
  fem.ls.bc[npIdP0] = 1
  fem.ls.bc[npIdP1] = 2
  fem.value[:] = 0.5
  fem.value[npIdP0] = 0.0
  fem.value[npIdP1] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  vis_color = dfm2.VisFEM_ColorContour(fem,name_color="value")
  vis_color.set_color_minmax()
  axis = dfm2.AxisXYZ(1.0)
  dfm2.glfw.winDraw3d([vis_color,axis])

def diffuse(msh,npIdP0,npIdP1):
  fem = dfm2.FEM_Diffuse(msh)
  fem.ls.bc[npIdP1] = 1
  fem.value[:] = 0.0
  fem.value[npIdP1] = 1.0
  ####
  vis_color = dfm2.VisFEM_ColorContour(fem,name_color="value")
  vis_color.draw_val_min = 0.0
  vis_color.draw_val_max = 1.0
  axis = dfm2.AxisXYZ(1.0)
  dfm2.glfw.winDraw3d([fem,vis_color,axis])


def linear_solid_static(msh,npIdP):
  fem = dfm2.FEM_SolidLinearStatic(msh,gravity=[0.3,0,0])
  fem.ls.bc[npIdP,:] = 1
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  vis_disp = dfm2.VisFEM_ColorContour(fem,name_disp="vec_val")
  axis = dfm2.AxisXYZ(1.0)
  dfm2.glfw.winDraw3d([vis_disp,axis])
  vis_disp.write_vtk("linearsolid3d.vtk")


def linear_solid_dynamic(msh,npIdP):
  fem = dfm2.FEM_SolidLinearDynamic(msh,gravity=[0.3,0,0])
  fem.ls.bc[npIdP,:] = 1
  ####
  vis_disp = dfm2.VisFEM_ColorContour(fem,name_disp="vec_val")
  axis = dfm2.AxisXYZ(1.0)
  dfm2.glfw.winDraw3d([fem,vis_disp,axis])


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