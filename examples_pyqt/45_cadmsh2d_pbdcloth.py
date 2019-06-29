####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys, math

from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtWidgets import QVBoxLayout

sys.path.append("..")
import pydelfem2 as dfm2
import pydelfem2.gl
import pydelfem2.pyqt
import pydelfem2.cadmshsim

from PyQt5.QtCore import Qt, QSize
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout


class Window_Pbd2D(dfm2.pyqt.QW_CadMshFem):
  def __init__(self):
    super(Window_Pbd2D, self).__init__()

    self.setWindowTitle("CAD_Mesh_Pbd2d")

    self.cadmsh = dfm2.cadmshsim.CadMesh2D_PBDCloth(edge_length=0.1)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, +0.8,+1, -0.8,+1, -1, +1])
    self.cadmsh.pbd.param_gravity_y = -0.1
    self.cadmsh.pbd.dt = 0.1

    self.ui_fem = dfm2.pyqt.QW_PBDCloth(self.cadmsh.pbd)

    super().init_UI()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window_Pbd2D()
  window.show()
  sys.exit(app.exec_())