####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys, math

from PyQt5.QtCore import Qt, QSize, pyqtSignal
from PyQt5.QtWidgets import QApplication, QWidget, QSizePolicy
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QCheckBox, QButtonGroup, QRadioButton

sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.gl
import delfem2.pyqt
import delfem2.cadmshfem


class Window_Poisson(dfm2.pyqt.QW_CadMshFem):
  def __init__(self):
    super(Window_Poisson, self).__init__()

    self.cadmsh = dfm2.cadmshfem.CadMesh2D_Diffuse(edge_length=0.05)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])

    self.ui_fem = dfm2.pyqt.QW_FEMDiffuse(self.cadmsh.fem)

    self.setWindowTitle("CAD_Mesh_Diffuse")

    super().init_UI()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window_Poisson()
  window.show()
  sys.exit(app.exec_())