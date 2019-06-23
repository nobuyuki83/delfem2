####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys, math

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout

sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.gl
import delfem2.pyqt
import delfem2.cadmshsim


class Window_SolidLinearEigen(dfm2.pyqt.QW_CadMshFem):
  def __init__(self):
    super(Window_SolidLinearEigen, self).__init__()
    self.setWindowTitle("CAD_Mesh_SolidLinearEigen")

    self.cadmsh = dfm2.cadmshsim.CadMesh2D_FEMSolidLinearEigen(edge_length=0.05)
    self.cadmsh.add_polygon([-1, -0.2, +1, -0.2, +1, +0.2, -1, +0.2])

    self.ui_fem = dfm2.pyqt.QW_FEMSolidLinearEigen(self.cadmsh.fem)

    super().init_UI()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window_SolidLinearEigen()
  window.show()
  sys.exit(app.exec_())