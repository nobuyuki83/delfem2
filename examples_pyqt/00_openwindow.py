####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys
import PySide2
from PySide2.QtWidgets import QApplication, QWidget


class MyWindow(QWidget): 
  def __init__(self):
    super().__init__()
    self.title = 'simple'
    self.width = 500
    self.height = 400
    self.initUI()

  def initUI(self):
    self.setWindowTitle(self.title)
    self.setGeometry(0, 0, self.width, self.height)
    self.show()

  def keyPressEvent(self, event):
    if event.text() == 'q':
      self.close()


def main():
  app = QApplication(sys.argv)
  gui = MyWindow()
  sys.exit(app.exec_())


if __name__ == "__main__":
  main()