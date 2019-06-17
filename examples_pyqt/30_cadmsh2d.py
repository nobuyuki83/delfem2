####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys

from PyQt5.QtWidgets import (QApplication, QVBoxLayout,
                             QWidget, QPushButton)

sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.gl
import delfem2.pyqt

class Window(QWidget):
  def __init__(self):
    super(Window, self).__init__()

    self.cadmsh = dfm2.CadMesh2D(0.2)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])

    self.glWidget = dfm2.pyqt.QOpenGLWidget_Cad2D()
    self.glWidget.cadobj = self.cadmsh

    self.ui_meshres = dfm2.pyqt.QUI_MeshRes(self.cadmsh)
    self.ui_meshres.func_updated = self.glWidget.update

    mainLayout = QVBoxLayout()
    mainLayout.addWidget(self.glWidget)
    mainLayout.addLayout(self.ui_meshres.hl)
    self.setLayout(mainLayout)

    self.setWindowTitle("CAD")

  def button_clicked(self,btn:QPushButton):
    if btn == self.btn:
      self.cadmsh.remesh()
      self.glWidget.update()

  def keyPressEvent(self, event):
    if event.text() == 'q':
      self.close()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window()
  window.show()
  sys.exit(app.exec_())
