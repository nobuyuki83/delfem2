import sys, numpy
sys.path.append("../module_py")
import dfm2

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
  fem.ls.vec_bc[npIdP0] = 1
  fem.ls.vec_bc[npIdP1] = 2
  fem.vec_val[:] = 0.5
  fem.vec_val[npIdP0] = 0.0
  fem.vec_val[npIdP1] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(msh,val_color=fem.vec_val[:,0])
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])

def diffuse(msh,npIdP0,npIdP1):
  fem = dfm2.FEM_Diffuse(msh)
  fem.ls.vec_bc[npIdP0] = 1
  fem.ls.vec_bc[npIdP1] = 2
  fem.vec_val[:] = 0.5
  fem.vec_val[npIdP0] = 0.0
  fem.vec_val[npIdP1] = 1.0
  ####
  field = dfm2.Field(msh,val_color=fem.vec_val[:,0])
  field.draw_val_min = 0.0
  field.draw_val_max = 1.0
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def linear_solid_static(msh,npIdP):
  fem = dfm2.FEM_LinearSolidStatic(msh,gravity=[0.3,0,0])
  fem.ls.vec_bc[npIdP,:] = 1
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(msh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])

def linear_solid_dynamic(msh,npIdP):
  fem = dfm2.FEM_LinearSolidDynamic(msh,gravity=[0.3,0,0])
  fem.ls.vec_bc[npIdP,:] = 1
  ####
  field = dfm2.Field(msh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


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