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

loop = [[0,0, 1,0, 1,1, 0,1]]
out = dfm2.triangulation(loop,edge_length=0.05)
msh = dfm2.Mesh()
msh.meshtri2d(out[0],out[1])
axis = dfm2.AxisXYZ(1)
dfm2.glfw.winDraw3d([msh, axis])

