import sys
sys.path.append("../module_py")
import dfm2


def poisson(cad,mesh):
  fem = dfm2.FEM_Poisson2D(mesh)
  dfm2.cad_setBCFlagEdge(fem.ls.vec_bc,
                         mesh.np_pos, [0, 1, 2, 3], cad, [0], 1, 1.0e-10)
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh,val_color=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])


def diffuse(cad,mesh):
  fem = dfm2.FEM_Diffuse2D(mesh)
  dfm2.cad_setBCFlagEdge(fem.ls.vec_bc,
                         mesh.np_pos, [0, 1, 2, 3], cad, [0], 1, 1.0e-10)
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
  dfm2.cad_setBCFlagEdge(fem.ls.vec_bc,
                         mesh.np_pos, [3], cad, [0,1], 1, 1.0e-10)
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([field,axis])


def linear_solid_dynamic(cad,mesh):
  fem = dfm2.FEM_LinearSolidDynamic2D(mesh)
  dfm2.cad_setBCFlagEdge(fem.ls.vec_bc,
                         mesh.np_pos, [3], cad, [0,1], 1, 1.0e-10)
  fem.solve()
  print(fem.ls.conv_hist)
  ####
  field = dfm2.Field(mesh,val_disp=fem.vec_val)
  axis = dfm2.AxisXYZ(1.0)
  dfm2.winDraw3d([fem,field,axis])


def main():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1.0])
  mesh = dfm2.mesh_cad(cad,0.05)
  #  dfm2.winDraw3d([cad,mesh])
  poisson(cad,mesh)
  diffuse(cad,mesh)
  linear_solid_static(cad,mesh)
  linear_solid_dynamic(cad,mesh)


if __name__ == "__main__":
  main()
