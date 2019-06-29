####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys

from PyQt5.QtWidgets import (QApplication, QVBoxLayout,
                             QWidget, QPushButton)

sys.path.append("..")
import pydelfem2 as dfm2
import pydelfem2.gl
import pydelfem2.pyqt

class Window(QWidget):
  def __init__(self):
    super(Window, self).__init__()

    self.cad = dfm2.Cad2D()
    self.cad.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])

    self.glWidget = dfm2.pyqt.QGLW_Cad2D()
    self.glWidget.cadobj = self.cad

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
