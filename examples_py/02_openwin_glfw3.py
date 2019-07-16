####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import sys
sys.path.append("..")
import PyDelFEM2 as dfm2
import PyDelFEM2.gl._glfw

msh = dfm2.Mesh()
msh.read("../test_inputs/bunny_2k.ply")
msh.scale_xyz(0.03)
dfm2.gl._glfw.winDraw3d([msh],winsize=(400,300))