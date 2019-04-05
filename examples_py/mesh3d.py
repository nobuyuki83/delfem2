from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2


def example1():
  sdf = dfm2.SDF()
  sdf.list_sdf.append( dfm2.SDF_Sphere(0.6,[-0.5,0,0],True) )
  sdf.list_sdf.append( dfm2.SDF_Sphere(0.6,[+0.5,0,0],True) )
  np_xyz,np_tet = dfm2.isosurface(sdf.list_sdf)
  print(np_xyz.shape,np_tet.shape)
  msh = dfm2.Mesh(np_xyz,np_tet)
  dfm2.winDraw3d([sdf,msh],winsize=(400,300))

if __name__ == "__main__":
  example1()