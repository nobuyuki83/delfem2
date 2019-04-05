from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

def make_mesh():
  sdf = dfm2.SDF()
  sdf.list_sdf.append( dfm2.SDF_Sphere(0.55,[-0.5,0,0],True) )
  sdf.list_sdf.append( dfm2.SDF_Sphere(0.55,[+0.5,0,0],True) )
  np_xyz,np_tet = dfm2.isosurface(sdf.list_sdf)
  print(np_xyz.shape,np_tet.shape)
  msh = dfm2.Mesh(np_xyz,np_tet,dfm2.Tet)
  return msh

def poission(msh):
  fem = dfm2.FEM_Poisson(msh)
  fem.vec_val[:] = 0.5
  np = msh.np_pos.shape[0]
  for ip in range(np):
    if msh.np_pos[ip,0] > 1:
      fem.ls.vec_bc[ip,0] = 1
      fem.vec_val[ip] = 0
    elif msh.np_pos[ip,0] < -1:
      fem.ls.vec_bc[ip,0] = 2
      fem.vec_val[ip,0] = 1
  print(fem.ls.conv_hist)
  fem.solve()
  ####
  field = dfm2.Field(msh,val_color=fem.vec_val[:,0])
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])

def diffuse(msh):
  fem = dfm2.FEM_Diffuse(msh)
  fem.vec_val[:] = 0.5
  np = msh.np_pos.shape[0]
  for ip in range(np):
    if msh.np_pos[ip,0] > 1:
      fem.ls.vec_bc[ip,0] = 1
      fem.vec_val[ip] = 0
    elif msh.np_pos[ip,0] < -1:
      fem.ls.vec_bc[ip,0] = 2
      fem.vec_val[ip,0] = 1
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(msh,val_color=fem.vec_val[:,0])
  field.draw_val_min = 0.0
  field.draw_val_max = 1.0
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def linear_solid_static(msh):
  fem = dfm2.FEM_LinearSolidStatic(msh,gravity=[0.3,0,0])
  for ip in range(msh.np_pos.shape[0]):
    if msh.np_pos[ip,0] < -1:
      fem.ls.vec_bc[ip,:] = 1
      fem.vec_val[ip,:] = 0.0
  print(fem.ls.conv_hist)
  fem.solve()
  ####
  field = dfm2.Field(msh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])

def linear_solid_dynamic(msh):
  fem = dfm2.FEM_LinearSolidDynamic(msh,gravity=[0.3,0,0])
  for ip in range(msh.np_pos.shape[0]):
    if msh.np_pos[ip,0] < -1:
      fem.ls.vec_bc[ip,:] = 1
      fem.vec_val[ip,:] = 0.0
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(msh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])

if __name__ == "__main__":
  msh = make_mesh()
  poission(msh)
  diffuse(msh)
  linear_solid_static(msh)
  linear_solid_dynamic(msh)