####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import OpenGL.GL as gl
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def main():
  msh = dfm2.Mesh()
  msh.read("../test_inputs/bunny_2k.ply")
  msh.color_face = [1,1,1,1]

  aabb = dfm2.AABB3( msh.minmax_xyz() )
  axis = dfm2.gl.AxisXYZ(100.0)
  axis.line_width=3

  sampler = dfm2.gl.CppGPUSampler()
  sampler.set_texture_property(size_res_width=256,size_res_height=256, is_rgba_8ui=True)
  sampler.set_coordinate(len_grid=0.2, depth_max=100.0,
                       org=[0, -50, 0], dir_prj=[0, -1.0, 0], dir_width=[1, 0, 0])  
  sampler.point_color = [0,0,1,1]
  sampler.len_axis = 10

  with dfm2.gl.glfw.GPUSamplerBufferGLFW([512, 512], format_color="4byte", is_depth=True):
    sampler.init_gl()
    sampler.start()
    gl.glClearColor(1,0,1,1)
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    gl.glDisable(gl.GL_LIGHTING)
    gl.glDisable(gl.GL_BLEND)
    msh.draw()
    sampler.end()
    sampler.get_depth()
    sampler.get_color()

#  np_depth = numpy.array(dfm2.depth_buffer(sampler),copy=True)
#  print(np_depth.shape)
#  numpy.savetxt("hoge.txt",np_depth)

  dfm2.gl.glfw.winDraw3d([msh,aabb,axis,sampler])
#  dfm2.gl.glfw.winDraw3d([msh,aabb,axis])

if __name__ == "__main__":
  main()