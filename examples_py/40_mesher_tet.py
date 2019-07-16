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

def example1():
  sdf = dfm2.SDF()
  sdf.list_sdf.append( dfm2.CppSDF_Sphere(0.6,[-0.5,0,0],True) )
  sdf.list_sdf.append( dfm2.CppSDF_Sphere(0.6,[+0.5,0,0],True) )
  dfm2.gl._glfw.winDraw3d([sdf],winsize=(400,300))

  np_xyz,np_tet = dfm2.isosurface(sdf.list_sdf)
  print(np_xyz.shape,np_tet.shape)
  msh = dfm2.Mesh(np_xyz,np_tet,dfm2.TET)
  dfm2.gl._glfw.winDraw3d([sdf,msh],winsize=(400,300))

if __name__ == "__main__":
  example1()