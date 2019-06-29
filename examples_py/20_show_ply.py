####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import sys
sys.path.append("..")
import pydelfem2 as dfm2
import pydelfem2.gl._glfw


def main():
  msh = dfm2.Mesh()
  msh.read("../test_inputs/bunny_2k.ply")
  aabb = dfm2.AABB3( msh.minmax_xyz() )
  dfm2.gl._glfw.winDraw3d([msh,aabb])

  msh.read("../test_inputs/bunny_1k.obj")
  aabb = dfm2.AABB3( msh.minmax_xyz() )
  dfm2.gl._glfw.winDraw3d([msh,aabb])

if __name__ == "__main__":
  main()