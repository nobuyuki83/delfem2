import sys
sys.path.append("../module_py")
import dfm2


def poisson(cad,mesh):
  fem = dfm2.FEM_Poisson(mesh,source=1.0)
  npIdP = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.vec_bc[npIdP] = 1
  fem.solve()
  field = dfm2.Field(mesh,val_color=fem.vec_val[:,0])
  dfm2.winDraw3d([field])
  field.write_vtk("poisson.vtk")


def poisson2(cad,mesh):
  fem = dfm2.FEM_Poisson(mesh,source=1.0)
  npIdP_fix = cad.points_edge([2], mesh.np_pos)
  fem.ls.vec_bc[npIdP_fix] = 1
  npIdP_ms = cad.points_edge([0], mesh.np_pos)
  fem.ls.vec_ms[:] = -1
  fem.ls.vec_ms[npIdP_ms] = npIdP_ms[0]
  fem.ls.vec_ms[npIdP_ms[0]] = -1
  fem.solve()
  field = dfm2.Field(mesh,val_color=fem.vec_val[:,0])
  dfm2.winDraw3d([field])


def diffuse(cad,mesh):
  fem = dfm2.FEM_Diffuse(mesh,source=1.0)
  npIdP = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.vec_bc[npIdP] = 1
  ####
  field = dfm2.Field(mesh,val_color=fem.vec_val)
  field.draw_val_min = 0.0
  field.draw_val_max = 0.3
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def linear_solid_static(cad,mesh):
  fem = dfm2.FEM_LinearSolidStatic(mesh,gravity=[0,-0.1])
  npIdP = cad.points_edge([3], mesh.np_pos)
  fem.ls.vec_bc[npIdP,:] = 1
  fem.solve()
  field = dfm2.Field(mesh,val_disp=fem.vec_val)
  dfm2.winDraw3d([field])


def linear_solid_dynamic(cad,mesh):
  fem = dfm2.FEM_LinearSolidDynamic(mesh,gravity=[0,-0.1])
  npIdP = cad.points_edge([3], mesh.np_pos)
  fem.ls.vec_bc[npIdP,:] = 1
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def storks_static(cad,mesh):
  fem = dfm2.FEM_StorksStatic2D(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.vec_bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh, val_color=fem.vec_val[:,2], val_disp=fem.vec_val[:,:2], disp_mode='hedgehog')
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])


def storks_dynamic(cad,mesh):
  fem = dfm2.FEM_StorksDynamic2D(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.vec_bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh, val_color=fem.vec_val[:,2], val_disp=fem.vec_val[:,:2], disp_mode='hedgehog')
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def navir_storks(cad,mesh):
  fem = dfm2.FEM_NavierStorks2D(mesh)
  npIdP0 = cad.points_edge([0,1,2,3], mesh.np_pos)
  fem.ls.vec_bc[npIdP0,0:2] = 1
  npIdP1 = cad.points_edge([2], mesh.np_pos)
  fem.vec_val[npIdP1,0] = 1.0
  ####
  field = dfm2.Field(mesh, val_color=fem.vec_val[:,2], val_disp=fem.vec_val[:,:2], disp_mode='hedgehog')
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def cloth(cad,mesh):
  fem = dfm2.FEM_Cloth(mesh)
  npIdP = cad.points_edge([2], mesh.np_pos)
  fem.ls.vec_bc[npIdP,0:3] = 1
  ####
  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=mesh.np_elm)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,mesh2,axis])


def main():
  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = cad.mesh(0.05)
  #  dfm2.winDraw3d([cad,mesh])
  poisson(cad,mesh)
  diffuse(cad,mesh)
  linear_solid_static(cad,mesh)
  linear_solid_dynamic(cad,mesh)
  storks_static(cad,mesh)
  storks_dynamic(cad,mesh)
  navir_storks(cad,mesh)
  cloth(cad,mesh)

  cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,0, +0,+0, 0,+1, -1,+1.0])
  mesh = cad.mesh(0.05)
  poisson2(cad,mesh)

if __name__ == "__main__":
  main()
