import sys
sys.path.append("../module_py")
import dfm2


def poisson(cad,mesh):
  fem = dfm2.FEM_Poisson2D(mesh)
  npIdP = dfm2.cad_getPointsEdge(cad,[0,1,2,3], mesh.np_pos, 1.0e-10)
  fem.ls.vec_bc[npIdP] = 1
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh,val_color=fem.vec_val[:,0])
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])


def diffuse(cad,mesh):
  fem = dfm2.FEM_Diffuse2D(mesh)
  npIdP = dfm2.cad_getPointsEdge(cad,[0,1,2,3], mesh.np_pos, 1.0e-10);
  fem.ls.vec_bc[npIdP] = 1
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh,val_color=fem.vec_val)
  field.draw_val_min = 0.0
  field.draw_val_max = 0.3
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def linear_solid_static(cad,mesh):
  fem = dfm2.FEM_LinearSolidStatic2D(mesh)
  npIdP = dfm2.cad_getPointsEdge(cad,[3], mesh.np_pos, 1.0e-10)
  fem.ls.vec_bc[npIdP,:] = 1
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])


def linear_solid_dynamic(cad,mesh):
  fem = dfm2.FEM_LinearSolidDynamic2D(mesh)
  npIdP = dfm2.cad_getPointsEdge(cad,[3], mesh.np_pos, 1.0e-10)
  fem.ls.vec_bc[npIdP,:] = 1
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def storks_static(cad,mesh):
  fem = dfm2.FEM_StorksStatic2D(mesh)
  npIdP0 = dfm2.cad_getPointsEdge(cad,[0,1,2,3], mesh.np_pos, 1.0e-10)
  fem.ls.vec_bc[npIdP0,0:2] = 1
  npIdP1 = dfm2.cad_getPointsEdge(cad,[2], mesh.np_pos, 1.0e-10)
  fem.vec_val[npIdP1,0] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh, val_color=fem.vec_val[:,2], val_disp=fem.vec_val[:,:2], disp_mode='hedgehog')
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])


def storks_dynamic(cad,mesh):
  fem = dfm2.FEM_StorksDynamic2D(mesh)
  npIdP0 = dfm2.cad_getPointsEdge(cad,[0,1,2,3], mesh.np_pos, 1.0e-10)
  fem.ls.vec_bc[npIdP0,0:2] = 1
  npIdP1 = dfm2.cad_getPointsEdge(cad,[2], mesh.np_pos, 1.0e-10)
  fem.vec_val[npIdP1,0] = 1.0
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh, val_color=fem.vec_val[:,2], val_disp=fem.vec_val[:,:2], disp_mode='hedgehog')
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def navir_storks(cad,mesh):
  fem = dfm2.FEM_NavierStorks2D(mesh)
  npIdP0 = dfm2.cad_getPointsEdge(cad,[0,1,2,3], mesh.np_pos, 1.0e-10)
  fem.ls.vec_bc[npIdP0,0:2] = 1
  npIdP1 = dfm2.cad_getPointsEdge(cad,[2], mesh.np_pos, 1.0e-10)
  fem.vec_val[npIdP1,0] = 1.0
  ####
  field = dfm2.Field(mesh, val_color=fem.vec_val[:,2], val_disp=fem.vec_val[:,:2], disp_mode='hedgehog')
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def cloth(cad,mesh):
  fem = dfm2.FEM_Cloth(mesh)
  npIdP = dfm2.cad_getPointsEdge(cad,[2], mesh.np_pos, 1.0e-10)
  fem.ls.vec_bc[npIdP,0:3] = 1
  ####
  print(fem.vec_val.shape)
  mesh2 = dfm2.Mesh(np_pos=fem.vec_val,np_elm=mesh.np_elm)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,mesh2,axis])


def main():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = dfm2.mesh_cad(cad,0.1)
  #  dfm2.winDraw3d([cad,mesh])
  poisson(cad,mesh)
  diffuse(cad,mesh)
  linear_solid_static(cad,mesh)
  linear_solid_dynamic(cad,mesh)
  storks_static(cad,mesh)
  storks_dynamic(cad,mesh)
  navir_storks(cad,mesh)
  cloth(cad,mesh)

if __name__ == "__main__":
  main()
