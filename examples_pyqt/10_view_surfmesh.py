####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys

from PySide2.QtCore import QSize
from PySide2.QtWidgets import (QApplication, QVBoxLayout, QOpenGLWidget,
                             QWidget, QPushButton)

import OpenGL.GL as gl

sys.path.append("..")
import PyDelFEM2 as dfm2
import PyDelFEM2.gl
import PyDelFEM2.qt

class Window(QWidget):
  def __init__(self):
    super(Window, self).__init__()

    self.msh = dfm2.Mesh()
    self.msh.read("../test_inputs/bunny_2k.ply");
    self.msh.scale_xyz(0.03)

    self.glWidget = GLWidget()
    self.glWidget.msh = self.msh

    self.btn = QPushButton('Button', self)

    mainLayout = QVBoxLayout()
    mainLayout.addWidget(self.glWidget)
    mainLayout.addWidget(self.btn)
    self.setLayout(mainLayout)

    self.setWindowTitle("CAD")

  def keyPressEvent(self, event):
    if event.text() == 'q':
      self.close()

class GLWidget(QOpenGLWidget):
  def __init__(self, parent=None):
    super(GLWidget, self).__init__(parent)
    self.msh = None
    self.nav = dfm2.qt.NavigationPyQt(view_height=1.0)

  def minimumSizeHint(self):
    return QSize(300, 300)

  def sizeHint(self):
    return QSize(600, 600)

  def initializeGL(self):
    print(dfm2.gl.getOpenglInfo())
    gl.glClearColor(0.8, 0.8, 1.0, 1.0)
    gl.glShadeModel(gl.GL_FLAT)
    gl.glEnable(gl.GL_DEPTH_TEST)
    gl.glEnable(gl.GL_CULL_FACE)
    gl.glDisable(gl.GL_LIGHTING)
    gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
    gl.glPolygonOffset(1.1, 4.0 )

  def paintGL(self):
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    self.nav.camera.set_gl_camera()
    self.msh.draw()

  def resizeGL(self, width, height):
    gl.glViewport(0,0,width,height)

  def mousePressEvent(self, event):
    self.nav.mouse(event, self.frameSize())

  def mouseMoveEvent(self, event):
    self.nav.motion(event, self.frameSize())
    self.update()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window()
  window.show()
  sys.exit(app.exec_())
