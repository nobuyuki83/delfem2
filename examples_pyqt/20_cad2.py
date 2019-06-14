####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys

from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import (QApplication, QVBoxLayout, QOpenGLWidget, QMenu,
                             QWidget, QPushButton)

import OpenGL.GL as gl

sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.gl
import delfem2.pyqt

class Window(QWidget):
  def __init__(self):
    super(Window, self).__init__()

    self.cad = dfm2.Cad2D()
    self.cad.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])

    self.glWidget = GLWidget()
    self.glWidget.cad = self.cad

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
    self.cad = None
    self.nav = dfm2.pyqt.NavigationPyQt(view_height=2.0)

  def minimumSizeHint(self):
    return QSize(300, 300)

  def sizeHint(self):
    return QSize(600, 600)

  def initializeGL(self):
    print(delfem2.gl.getOpenglInfo())
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
    self.cad.draw()

  def resizeGL(self, width, height):
    gl.glViewport(0,0,width,height)

  def mousePressEvent(self, event):
    self.makeCurrent()
    self.nav.mouse(event, self.frameSize())
    src,dir = self.nav.mouse_src_dir()
    self.cad.cad.pick(src[0],src[1],self.nav.camera.view_height)
    self.update()

  def mouseReleaseEvent(self, event):
    self.cad.cad.ivtx_picked = -1
    self.update()

  def mouseMoveEvent(self, event):
    self.makeCurrent()
    self.nav.motion(event, self.frameSize())
    src0,src1,dir = self.nav.motion_src_dir()
    self.cad.motion(src0,src1,dir)
    self.update()

  def contextMenuEvent(self, event):
    src,dir = self.nav.mouse_src_dir()
    ###
    menu = QMenu(self)
    actionAddPoint = None
    if self.cad.cad.iedge_picked != -1:
      actionAddPoint = menu.addAction("add point")
    action = menu.exec_(self.mapToGlobal(event.pos()))
    ####
    if action == actionAddPoint != None:
      self.cad.add_point_edge(src[0],src[1],self.cad.cad.iedge_picked)
      self.update()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window()
  window.show()
  sys.exit(app.exec_())
