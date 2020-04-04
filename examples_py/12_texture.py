####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import sys, cv2, math
sys.path.append("..")
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def main():
  path_img = "../test_inputs/lenna.png"
  np_img = cv2.imread(path_img)
  axis = dfm2.gl.AxisXYZ(100)
  tex = dfm2.gl.get_texture(np_img,"rgb")
  dfm2.gl.glfw.winDraw3d([axis,tex],winsize=(400,300),
                          camera_rotation=[math.pi, 0, 0])

if __name__ == "__main__":
  main()