####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def main():
  cad = dfm2.Cad2D()
  cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  dfm2.gl.glfw.winDraw3d([cad])

  cad.add_vtx_face(0, [0.0, 0.0])
  dfm2.gl.glfw.winDraw3d([cad])

  cad.add_vtx_edge(2, [0.0, 0.8])
  dfm2.gl.glfw.winDraw3d([cad])

  cad.set_edge_type(0, dfm2.CAD_EDGE_GEOM_BEZIER_CUBIC, [0.2, 0.3, 0.8, 0.3])
  dfm2.gl.glfw.winDraw3d([cad])

  cad.set_edge_type(0, dfm2.CAD_EDGE_GEOM_LINE, [])
  dfm2.gl.glfw.winDraw3d([cad])

  cad.clear()
  cad.add_polygon(list_xy=[-1000,-1000, +1000,-1000, +1000,+1000, -1000,+1000])
  dfm2.gl.glfw.winDraw3d([cad])

  for itr in range(3):
    cad.clear()
    cad.import_svg("../test_inputs/shape"+str(itr)+".svg")
    dfm2.gl.glfw.winDraw3d([cad])


if __name__ == "__main__":
  main()
