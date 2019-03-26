import numpy
import sys
sys.path.append("../module_py")
import dfm2

def main():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])

  mesh = dfm2.mesh_cad(cad,0.05)

#  fem = dfm2.FEM(mesh)
  # matrix
  mat = dfm2.MatrixSquareSparse()
  mat.initialize(mesh.np_pos.shape[0],1, True)
  psup_ind, psup = mesh.psup()
  dfm2.sortIndexedArray(psup_ind, psup)
  dfm2.matrixSquareSparse_setPattern(mat,psup_ind,psup)

  # preconditioner
  mat_prec = dfm2.PreconditionerILU()
  dfm2.precond_ilu0(mat_prec,mat)

  # vectors
  np = mesh.np_pos.shape[0]
  print("number of points",np)
  vec_val = numpy.zeros((np,),dtype=numpy.float64) # initial guess is zero
  vec_bc = numpy.zeros((np,),dtype=numpy.int32)
  vec_bc[0] = 1
  dfm2.cad_setBCFlagEdge(vec_bc,mesh.np_pos,[0,1,2,3],cad,1,1.0e-10)

  ####
  mat.setZero()
  vec_f = numpy.ndarray((np,), dtype=numpy.float64)
  vec_f[:] = 0.0
  dfm2.mergeLinSys_poission2D(mat,vec_f,
                              1.0,0.1,
                              mesh.np_pos,mesh.np_elm,vec_val)
  print(numpy.linalg.norm(vec_f))
  #### setting bc
  vec_f[vec_bc!=0] = 0.0
  dfm2.matrixSquareSparse_setFixBC(mat,vec_bc)

  #### solving matrix
  mat_prec.set_value(mat)
  mat_prec.ilu_decomp()
  vec_x = numpy.zeros((np,), dtype=numpy.float64)
  dfm2.linsys_solve_pcg(vec_f,vec_x,
                        0.0001, 100, mat, mat_prec)
  vec_x[vec_bc!=0] = 0.0
  vec_val += vec_x
  print(vec_val.min(),vec_val.max())

  field = dfm2.Field(mesh,vec_val)

  axis = dfm2.AxisXYZ(1.0)

#  dfm2.winDraw3d([cad,mesh])
  dfm2.winDraw3d([field,axis])



if __name__ == "__main__":
  main()
