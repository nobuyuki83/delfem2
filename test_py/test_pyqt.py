####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import pytest
import os, sys
from PySide2.QtCore import QTimer
from PySide2.QtWidgets import QApplication, QWidget, QVBoxLayout

import PyDelFEM2 as dfm2
import PyDelFEM2.qt
os.chdir(os.path.dirname(os.path.abspath(__file__))) # for python3 setup.py test

app = QApplication(sys.argv)


@pytest.mark.qttest
class Test_OpenWindow():
  class MyWindow(QWidget):
    def __init__(self):
      super().__init__()
      self.setWindowTitle("simple window test")
      self.setGeometry(0, 0, 500, 400)

  timer = QTimer()
  timer.timeout.connect(lambda: app.closeAllWindows())

  def test0(self):
    global app
    self.timer.start(1000) # 1 sec
    gui = self.MyWindow()
    gui.show()
    app.exec_()

'''
class Test_QGLW_Mesh(unittest.TestCase):
  class Window(QWidget):
    def __init__(self):
      super().__init__()
      self.msh = dfm2.Mesh()
      self.msh.read("../test_inputs/bunny_2k.ply");
      self.msh.scale_xyz(0.03)
      self.glWidget = dfm2.qt.QGLW_Mesh()
      self.glWidget.msh = self.msh
      mainLayout = QVBoxLayout()
      mainLayout.addWidget(self.glWidget)
      self.setLayout(mainLayout)
      self.setWindowTitle("CAD")

  timer = QTimer()
  timer.timeout.connect(lambda: app.closeAllWindows())

  def test0(self):
    global app
    gui = dfm2.qt.QGLW_Mesh()
    gui.show()
    self.timer.start(1000) # 1 sec
    app.exec_()

  def test1(self):
    global appn
    print("hoge0")
    gui = dfm2.qt.QGLW_Mesh()
    print("hoge1")
    gui.msh = dfm2.Mesh()
    print("hoge2")
    gui.msh.read("../test_inputs/bunny_2k.ply")
    print("hoge3")
    gui.msh.scale_xyz(0.03)
    print("hoge4")
    gui.show()
    print("hoge5")
    self.timer.start(1000) # 1 sec
    print("hoge6")
    app.exec_()

  def test2(self):
    global app
    gui = self.Window()
    gui.show()
    self.timer.start(1000) # 1 sec
    app.exec_()


class Test_QGLW_CAD2D(unittest.TestCase):
  class QW_Cad2D(QWidget):
    def __init__(self):
      super().__init__()

      self.cad = dfm2.Cad2D()
      self.cad.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])

      self.glWidget = dfm2.qt.QGLW_Cad2D()
      self.glWidget.cadobj = self.cad

      mainLayout = QVBoxLayout()
      mainLayout.addWidget(self.glWidget)
      self.setLayout(mainLayout)

      self.setWindowTitle("CAD")

  timer = QTimer()
  timer.timeout.connect(lambda: app.closeAllWindows())

  def test0(self):
    global app
    self.timer.start(1000) # 1 sec
    gui = dfm2.qt.QGLW_Cad2D()
    gui.cadobj = dfm2.Cad2D()
    gui.cadobj.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])
    gui.show()
    app.exec_()

  def test1(self):
    global app
    self.timer.start(1000) # 1 sec
    gui = self.QW_Cad2D()
    gui.show()
    app.exec_()

'''    