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


def example0():
  dmsh = dfm2.MeshDynTri2D()
  dmsh.meshing_loops([[0,0, 1,0, 1,1, 0,1]], edge_length=0.05)
  axis = dfm2.AxisXYZ(1)
  dfm2.glfw.winDraw3d([dmsh, axis])
  #####
  msh = dfm2.Mesh()
  msh.set_extrude(dmsh,1)
  msh.np_pos[:,2] *= 0.1
  dfm2.glfw.winDraw3d([msh, axis])


if __name__ == "__main__":
  example0()


