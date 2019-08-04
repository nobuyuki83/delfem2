####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def example1():
  sdf0 = dfm2.CppSDF3_Sphere(0.6,[-0.5,0,0],True)
  sdf1 = dfm2.CppSDF3_Sphere(0.6,[+0.5,0,0],True)
  dfm2.gl.glfw.winDraw3d([sdf0,sdf1],winsize=(400,300))

  np_xyz,np_tet = dfm2.isosurface([sdf0,sdf1])
  print(np_xyz.shape,np_tet.shape)
  msh = dfm2.Mesh(np_xyz,np_tet,dfm2.TET)
  dfm2.gl.glfw.winDraw3d([sdf0,sdf1,msh],winsize=(400,300))

if __name__ == "__main__":
  example1()