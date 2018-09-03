from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2


def main():
  mshelm = dfm2.MeshElem("../test_inputs/bunny_2k.ply");
  aabb = dfm2.AABB3( mshelm.minmax_xyz() )
  axis = dfm2.AxisXYZ(100.0)
  axis.line_width=3

  depth = dfm2.Depth(size_res_width=256,size_res_height=256,len_grid=0.5,depth_max=100.0,
                      org=[0,-50,0],dir_prj=[0,-1.0,0],dir_width=[1,0,0])
  depth.color = [0,0,1]

  depth_context = dfm2.DepthContext(win_size=[512,512])
  dfm2.take_depth_shot(mshelm.draw,depth,depth_context)
  depth_context.close()

  dfm2.winDraw3d([mshelm,aabb,depth,axis])

if __name__ == "__main__":
  main()