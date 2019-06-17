####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import sys
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw
import delfem2.cadmshfem

def mesh():
  cmf = dfm2.CadMesh2D(edge_length=0.1)
  cmf.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  dfm2.glfw.winDraw3d([cmf])

def poisson():
  cmf = dfm2.cadmshfem.CadMesh2D_Poisson(edge_length=0.1)
  cmf.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  dfm2.glfw.winDraw3d([cmf])

if __name__ == "__main__":
  mesh()
  poisson()