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
timer = QTimer()
timer.timeout.connect(lambda: app.closeAllWindows())

class Test_OpenWindow():
  class MyWindow(QWidget):
    def __init__(self):
      super().__init__()
      self.setWindowTitle("simple window test")
      self.setGeometry(0, 0, 500, 400)

  def test0(self):
    global app, timer
    timer.start(1000) # 1 sec
    gui = self.MyWindow()
    gui.show()
    app.exec_()


class Test_QGLW_Mesh():
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

  def test0(self):
    global app, timer
    timer.start(1000) # 1 sec    
    gui = dfm2.qt.QGLW_Mesh()
    gui.show()
    app.exec_()

  def test1(self):
    global app, timer
    gui = dfm2.qt.QGLW_Mesh()    
    gui.msh = dfm2.Mesh()
    gui.msh.read("../test_inputs/bunny_2k.ply")
    gui.msh.scale_xyz(0.03)    
    timer.start(1000) # 1 sec    
    gui.show() # segment fault here
    app.exec_()

  def test2(self):
    global app, timer
    gui = self.Window()
    timer.start(1000) # 1 sec        
    gui.show()
    app.exec_()


class Test_QGLW_CAD2D():
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

  def test0(self):
    global app, timer
    gui = dfm2.qt.QGLW_Cad2D()
    gui.cadobj = dfm2.Cad2D()
    gui.cadobj.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])
    timer.start(1000) # 1 sec             
    gui.show()
    app.exec_()

  def test1(self):
    global app
    gui = self.QW_Cad2D()
    timer.start(1000) # 1 sec                
    gui.show()
    app.exec_()
