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

import numpy, math

def main():
  A0 = numpy.zeros((16,32))
  for iy in range(A0.shape[0]):
    for ix in range(A0.shape[1]):
      A0[iy,ix] = 5*math.sin(ix*0.4)*math.cos(iy*0.6) + 5

  msh = dfm2.mesh_grid((A0.shape[1],A0.shape[0]))
  msh.np_pos = numpy.hstack((msh.np_pos,A0.reshape((-1,1))))

  axis = dfm2.AxisXYZ(32)
  dfm2.glfw.winDraw3d([msh,axis],(400,400))


if __name__ == "__main__":
  main()