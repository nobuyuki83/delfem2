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


class Window_Poisson(QWidget):
  def __init__(self):
    super(Window_Poisson, self).__init__()
    self.setWindowTitle("CAD_Mesh_Poisson")

    self.cadmsh = dfm2.cadmshfem.CadMesh2D_Poisson(edge_length=0.05)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])

    self.init_UI()

  def sizeHint(self):
    return QSize(1200, 600)

  def init_UI(self):
    self.glWidget0 = dfm2.pyqt.QOpenGLWidget_Cad2D()
    self.glWidget0.cadobj = self.cadmsh
    self.glWidget0.updated = self.updated

    self.glWidget1 = dfm2.pyqt.QOpenGLWidget_Cad2D()
    self.glWidget1.cadobj = self.cadmsh
    self.glWidget1.updated = self.updated

    self.ui_meshres = dfm2.pyqt.QUI_MeshRes(self.cadmsh)
    self.ui_meshres.func_updated = self.updated

    self.ui_poisson = dfm2.pyqt.QUI_FEMPoisson(self.cadmsh.fem)
    self.ui_poisson.func_updated = self.updated

    self.ui_numwin = dfm2.pyqt.QUI_NumWin()
    self.ui_numwin.changed.connect(self.numwin_changed)

    self.layout_param = QVBoxLayout()
    self.layout_param.addWidget(self.ui_numwin,  alignment=Qt.AlignLeft)
    self.layout_param.addWidget(self.ui_meshres, alignment=Qt.AlignLeft)
    self.layout_param.addWidget(self.ui_poisson, alignment=Qt.AlignLeft)
    self.layout_param.addStretch()

    self.mainLayout = QHBoxLayout()
    self.mainLayout.addWidget(self.glWidget0, alignment=Qt.AlignLeft, stretch=2)
    self.mainLayout.addWidget(self.glWidget1, alignment=Qt.AlignLeft, stretch=2)
    self.mainLayout.addLayout(self.layout_param, stretch=0)
    self.setLayout(self.mainLayout)

  def updated(self):
    self.glWidget0.update()
    self.glWidget1.update()

  def keyPressEvent(self, event):
    if event.text() == 'q':
      self.close()

  def numwin_changed(self,numwin:int):
    if numwin == 1:
      self.glWidget1.hide()
    elif numwin == 2:
      self.glWidget1.show()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window_Poisson()
  window.show()
  sys.exit(app.exec_())