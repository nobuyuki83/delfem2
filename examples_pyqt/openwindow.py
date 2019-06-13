import sys
import PyQt5
from PyQt5.QtWidgets import *


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


def main():
  app = QApplication(sys.argv)
  gui = MyWindow()
  sys.exit(app.exec_())


if __name__ == "__main__":
  main()