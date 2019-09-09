####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import os, sys
import pytest
from time import sleep

from PySide2.QtCore import Qt, QPointF, QEvent
from PySide2.QtCore import QTimer
from PySide2.QtWidgets import QApplication, QWidget, QVBoxLayout
from PySide2.QtGui import QMouseEvent

import PyDelFEM2 as dfm2
import PyDelFEM2.qt

os.chdir(os.path.dirname(os.path.abspath(__file__))) # for python3 setup.py test

app = QApplication(sys.argv)

def test_OpenWindow():
  class MyWindow(QWidget):
    def __init__(self):
      super().__init__()
      self.setWindowTitle("simple window test")
      self.setGeometry(0, 0, 500, 400)

  gui = MyWindow()
  icnt = 0
  def events():
    nonlocal icnt, gui
    icnt = icnt+1
    if icnt >= 30:
      gui.close()

  timer = QTimer()
  timer.setInterval(30)
  timer.timeout.connect(events)
  timer.start()
  gui.show()
  app.exec_() #event loop


def test_QGLWMesh():
  gui = dfm2.qt.QGLW_Mesh()
  icnt = 0

  def events():
    nonlocal icnt, gui
    icnt = icnt+1
    if icnt == 10:
      gui.msh = dfm2.Mesh()
      gui.msh.read("../test_inputs/bunny_2k.ply")
      gui.msh.scale_xyz(0.03)
    if icnt == 20:
      gui.make_buffer()
      event = QMouseEvent(QEvent.MouseButtonPress, QPointF(200.0,200.0), Qt.LeftButton, Qt.NoButton, Qt.ShiftModifier)
      gui.mousePressEvent(event)
    if icnt > 20 and icnt < 30:
      event = QMouseEvent(QEvent.MouseMove, QPointF(200.0+(icnt-20)*5,200.0), Qt.LeftButton, Qt.NoButton, Qt.NoModifier)
      gui.mouseMoveEvent(event)
    if icnt >= 30:
      gui.close()
    gui.update()

  gui.show()
  timer = QTimer()
  timer.setInterval(30)
  timer.timeout.connect(events)
  timer.start() # 1 sec
  app.exec_() # event loop



def test_QGLW_Cad2D():
  gui = dfm2.qt.QGLW_Cad2D()
  icnt = 0

  def events():
    nonlocal icnt, gui
    if icnt == 10:
      gui.cadobj = dfm2.Cad2D()
      gui.cadobj.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])
    if icnt == 20:
      gui.cadobj.pick(-1,+1, view_height=2)
      assert gui.cadobj.ivtx_picked() == 3
    if icnt > 20 and icnt < 30:
      gui.cadobj.motion([-1,+1,0],
                        [-1 - (icnt - 20) * 0.02, +1, 0], [0, 0, 1])
    if icnt == 30:
      gui.close()
    icnt += 1
    gui.update()

  gui.show()
  timer = QTimer()
  timer.setInterval(50)
  timer.timeout.connect(events)
  timer.start() # 1 sec
  app.exec_() # loop