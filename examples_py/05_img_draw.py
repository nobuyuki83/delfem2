####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import matplotlib.pyplot as plt

import PyDelFEM2 as dfm2
import PyDelFEM2.gl._glfw

msh = dfm2.Mesh()
msh.read("../test_inputs/bunny_2k.ply")
img = dfm2.gl._glfw.imgDraw3d([msh],winsize=(400,300))
plt.imshow(img)
plt.show()
