####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import PyDelFEM2 as dfm2
import PyDelFEM2.gl._glfw

def example0():
  dmsh = dfm2.MeshDynTri2D()
  dmsh.meshing_loops([[0,0, 1,0, 1,1, 0,1]], edge_length=0.05)
  axis = dfm2.gl.AxisXYZ(1)
  dfm2.gl._glfw.winDraw3d([dmsh, axis], camera_orientation=[+1,+1,-1, 0,1,0])
  #####
  msh = dfm2.Mesh()
  msh.set_extrude(dmsh,10)
  msh.np_pos[:,2] *= 0.1
  dfm2.gl._glfw.winDraw3d([msh, axis], camera_orientation=[+1,+1,-1, 0,1,0])


if __name__ == "__main__":
  example0()


