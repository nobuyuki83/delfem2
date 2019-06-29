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


class Window_SolidLinearStatic(dfm2.pyqt.QW_CadMshFem):
  def __init__(self):
    super(Window_SolidLinearStatic, self).__init__()

    self.setWindowTitle("CAD_Mesh_SolidLinearStatic")

    self.cadmsh = dfm2.cadmshsim.CadMesh2D_FEMSolidLinearStatic(edge_length=0.05)
    self.cadmsh.fem.param_gravity_y = -0.1
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])

    self.ui_fem = dfm2.pyqt.QW_FEMSolidLinearStatic(self.cadmsh.fem)

    super().init_UI()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window_SolidLinearStatic()
  window.show()
  sys.exit(app.exec_())