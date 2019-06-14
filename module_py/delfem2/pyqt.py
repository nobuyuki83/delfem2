####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys

from PyQt5.QtCore import pyqtSignal, QPoint, QSize, Qt, QEvent
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QApplication, QHBoxLayout, QOpenGLWidget, QSlider,
                             QWidget, QPushButton)

import OpenGL.GL as gl

sys.path.append("../module_py")
import numpy
import delfem2 as dfm2
import delfem2.gl

class NavigationPyQt:
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
    w0,h0 = size.width(),size.height()
    self.mouse_pre_x,self.mouse_pre_y = self.mouse_x, self.mouse_y
    self.mouse_x = (2.0 * event.x() - w0) / w0
    self.mouse_y = (h0 - 2.0 * event.y()) / h0
    if self.button & Qt.LeftButton:
      if self.modifiers & Qt.AltModifier:
        self.camera.rotation(self.mouse_x, self.mouse_y,
                             self.mouse_pre_x, self.mouse_pre_y)
      elif self.modifiers & Qt.ShiftModifier:
        self.camera.translation(self.mouse_x, self.mouse_y,
                                self.mouse_pre_x, self.mouse_pre_y)

  def btn(self):
    if self.button & Qt.LeftButton: return 0
    elif self.button & Qt.RightButton: return 1
    elif self.button & Qt.MiddleButton: return 2
    return -1

  def mouse_src_dir(self):
    self.camera.set_gl_camera()
    mMV = gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX)
    mPj = gl.glGetFloatv(gl.GL_PROJECTION_MATRIX)
    src = dfm2.gl.screenUnProjection(
      ([float(self.mouse_x), float(self.mouse_y), 0.0]),
      mMV, mPj)
    dir = dfm2.gl.screenUnProjectionDirection((0,0,1), mMV,mPj)
    return src,dir

  def motion_src_dir(self):
    self.camera.set_gl_camera()
    mMV = gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX)
    mPj = gl.glGetFloatv(gl.GL_PROJECTION_MATRIX)
    src0 = dfm2.gl.screenUnProjection(
      numpy.array([float(self.mouse_pre_x), float(self.mouse_pre_y), 0.0]),
      mMV, mPj)
    src1 = dfm2.gl.screenUnProjection(
      numpy.array([float(self.mouse_x), float(self.mouse_y), 0.0]),
      mMV, mPj)
    dir = dfm2.gl.screenUnProjectionDirection((0,0,1), mMV,mPj)
    return src0,src1,dir

def setClearColor(c:QColor):
  gl.glClearColor(c.redF(), c.greenF(), c.blueF(), c.alphaF())

def setColor(c:QColor):
  gl.glColor4f(c.redF(), c.greenF(), c.blueF(), c.alphaF())