####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys, math

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

    self.cadmsh = dfm2.CadMesh2D(0.2)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])

    self.glWidget = GLWidget()
    self.glWidget.cadobj = self.cadmsh

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
    self.cadobj = None
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
    self.cadobj.draw()

  def resizeGL(self, width, height):
    gl.glViewport(0,0,width,height)

  def mousePressEvent(self, event):
    self.makeCurrent()
    self.nav.mouse(event, self.frameSize())
    src,dir = self.nav.mouse_src_dir()
    self.cadobj.pick(src[0],src[1],self.nav.camera.view_height)
    self.update()

  def mouseReleaseEvent(self, event):
    self.cadobj.clean_picked()
    self.update()

  def mouseMoveEvent(self, event):
    self.makeCurrent()
    self.nav.motion(event, self.frameSize())
    src0,src1,dir = self.nav.motion_src_dir()
    self.cadobj.drag_picked(src1[0],src1[1], src0[0],src0[1])
    self.update()

  def contextMenuEvent(self, event):
    src,dir = self.nav.mouse_src_dir()
    ###
    menu = QMenu(self)
    actionDelVtx = None
    actionAddVtx = None
    actionAddSquare = None
    if self.cadobj.ivtx_picked() != -1:
      actionDelVtx = menu.addAction("delete vtx")
    elif self.cadobj.iedge_picked() != -1:
      actionAddVtx = menu.addAction("add vtx")
    else:
      actionAddSquare = menu.addAction("add square")
    action = menu.exec_(self.mapToGlobal(event.pos()))
    ####
    if action == actionAddVtx != None:
      self.cadobj.add_vtx_edge(src[0],src[1],self.cadobj.iedge_picked())
    if action == actionAddSquare != None:
      x0, y0 = src[0], src[1]
      self.cadobj.add_polygon([x0 - 1, y0 - 1, x0 + 1, y0 - 1, x0 + 1, y0 + 1, x0 - 1, y0 + 1])
    self.update()

  def wheelEvent(self,event):
    v = 0.1*(event.angleDelta().y()/120.0)
    self.nav.camera.scale *= math.exp(v)
    self.update()

if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = Window()
  window.show()
  sys.exit(app.exec_())
