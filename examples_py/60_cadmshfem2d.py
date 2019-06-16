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

def mesh():
  cmf = dfm2.CadMesh2D(edge_length=0.1)
  cmf.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
  dfm2.glfw.winDraw3d([cmf])


class CadMesh_Poisson(dfm2.CadMesh2D):

  def __init__(self,edge_length:float):
    super().__init__(edge_length)
    self.fem = dfm2.FEM_Poisson(self.msh,source=1.0)
    self.vis = dfm2.VisFEM_ColorContour(self.fem,name_color='value')
    self.remesh()

  def draw(self):
    self.ccad.draw()
    self.vis.draw()

  def minmax_xyz(self):
    return self.msh.minmax_xyz()

  def remesh(self):
    super().remesh()
    self.fem.updated_mesh()
    npIdP = dfm2.cad_getPointsEdge(self.ccad, [0, 1, 2, 3], self.msh.np_pos, 0.001)
    self.fem.value[npIdP] = 0.0
    self.fem.ls.bc[npIdP] = 1
    self.fem.solve()
    
  def motion(self,src0,src1,dir):
    super().motion(src0,src1,dir)
    self.fem.solve()


def poisson():
  cmf = CadMesh_Poisson(edge_length=0.1)
  cmf.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  dfm2.glfw.winDraw3d([cmf])

if __name__ == "__main__":
  mesh()
  poisson()