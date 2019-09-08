####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

from PySide2.QtCore import QSize
from PySide2.QtWidgets import (QApplication, QVBoxLayout, QOpenGLWidget,
                             QWidget, QPushButton)

import OpenGL.GL as gl

import PyDelFEM2 as dfm2
import PyDelFEM2.gl
import PyDelFEM2.qt

import sys

class Window(QWidget):
  def __init__(self):
    super(Window, self).__init__()

    self.msh = dfm2.Mesh()
    self.msh.read("../test_inputs/bunny_2k.ply");
    self.msh.scale_xyz(0.03)

    self.glWidget = dfm2.qt.QGLW_Mesh()
    self.glWidget.msh = self.msh

    mainLayout = QVBoxLayout()
    mainLayout.addWidget(self.glWidget)
    self.setLayout(mainLayout)

    self.setWindowTitle("CAD")

  def keyPressEvent(self, event):
    if event.text() == 'q':
      self.close()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window()
  window.show()
  sys.exit(app.exec_())
