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

def main():
  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  dfm2.gl._glfw.winDraw3d([cad])

  cad.add_vtx_edge(2, [0.0, 0.0])
  dfm2.gl._glfw.winDraw3d([cad])

  cad.set_edge_type(0, 1, [0.2, 0.3, -0.2, 0.3])
  dfm2.gl._glfw.winDraw3d([cad])

  cad.set_edge_type(0, 0, [])
  dfm2.gl._glfw.winDraw3d([cad])

if __name__ == "__main__":
  main()
