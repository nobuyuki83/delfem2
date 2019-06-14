import sys

from PyQt5.QtCore import pyqtSignal, QPoint, QSize, Qt, QEvent
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QApplication, QHBoxLayout, QOpenGLWidget, QSlider,
                             QWidget, QPushButton)

import OpenGL.GL as gl

sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.gl

class WindowManagerPyQt:
  def __init__(self,view_height):
    self.camera = dfm2.gl.Camera(view_height)
    self.mouse_x = 0.0
    self.mouse_y = 0.0
    self.mouse_pre_x = 0.0
    self.mouse_pre_y = 0.0
    self.button = Qt.MouseButton()
    self.modifiers = None

  def mouse(self,event,size:QSize):
    self.mouse_x = (2.0 * event.x() - size.width()) / size.width()
    self.mouse_y = (size.height() - 2.0 * event.y()) / size.height()
    self.button = event.button()
    self.modifiers = event.modifiers()

  def motion(self,event,size:QSize):
    self.mouse_pre_x,self.mouse_pre_y = self.mouse_x, self.mouse_y
    if self.button & Qt.LeftButton:
      if self.modifiers & Qt.AltModifier:
        self.mouse_x, self.mouse_y = self.camera.rotation(event.x(), event.y(),
                                                          self.mouse_x, self.mouse_y, size.width(), size.height())
      elif self.modifiers & Qt.ShiftModifier:
        self.mouse_x, self.mouse_y = self.camera.translation(event.x(), event.y(),
                                                          self.mouse_x, self.mouse_y, size.width(), size.height())
    self.mouse_x = (2.0 * event.x() - size.width()) / size.width()
    self.mouse_y = (size.height() - 2.0 * event.y()) / size.height()

def setClearColor(c:QColor):
  gl.glClearColor(c.redF(), c.greenF(), c.blueF(), c.alphaF())

def setColor(c:QColor):
  gl.glColor4f(c.redF(), c.greenF(), c.blueF(), c.alphaF())