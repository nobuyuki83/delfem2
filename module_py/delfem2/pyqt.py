####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys, math, numpy
import OpenGL.GL as gl

from PyQt5.QtCore import pyqtSignal, QPoint, QSize, Qt, QEvent
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QOpenGLWidget, QMenu, QWidget
from PyQt5.QtWidgets import QPushButton, QLabel, QSlider, QHBoxLayout, QCheckBox, QVBoxLayout, QLabel

sys.path.append("../module_py")
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


class QOpenGLWidget_Cad2D(QOpenGLWidget):
  def __init__(self, parent=None):
    super(QOpenGLWidget_Cad2D, self).__init__(parent)
    self.cadobj = None
    self.nav = dfm2.pyqt.NavigationPyQt(view_height=2.0)

  def minimumSizeHint(self):
    return QSize(300, 300)

  def sizeHint(self):
    return QSize(600, 600)

  def initializeGL(self):
    print(delfem2.gl.getOpenglInfo())
    gl.glClearColor(0.8, 0.8, 1.0, 1.0)
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
    self.cadobj.pick(src[0], src[1], self.nav.camera.view_height)
    self.update()

  def mouseReleaseEvent(self, event):
    self.cadobj.clean_picked()
    self.update()

  def mouseMoveEvent(self, event):
    self.makeCurrent()
    self.nav.motion(event, self.frameSize())
    src0,src1,dir = self.nav.motion_src_dir()
    self.cadobj.motion(src0, src1, dir)
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
      self.cadobj.add_vtx_edge(src[0], src[1], self.cadobj.ccad.iedge_picked)
    if action == actionAddSquare != None:
      x0,y0 = src[0],src[1]
      self.cadobj.add_polygon([x0-1, y0-1, x0+1, y0-1, x0+1, y0+1, x0-1, y0+1])
    self.update()

  def wheelEvent(self,event):
    v = 0.1*(event.angleDelta().y()/120.0)
    self.nav.camera.scale *= math.exp(v)
    self.update()


class QUI_MeshRes(QWidget):
  def __init__(self,cadmsh:delfem2.CadMesh2D):
    super(QUI_MeshRes, self).__init__()
    self.cadmsh = cadmsh
    self.func_updated = None

    self.btn = QPushButton('remesh', self)
    self.btn.clicked.connect(lambda: self.button_clicked(self.btn))

    self.lbl = QLabel("50", self)

    self.sp = QSlider(Qt.Horizontal)
    self.sp.setMinimum(10)
    self.sp.setMaximum(100)
    self.sp.setValue(50)
    self.sp.valueChanged.connect(self.slider_moved)

    self.b1 = QCheckBox("Sync Mesh To CAD")
    self.b1.setChecked(True)
    self.b1.stateChanged.connect(lambda: self.btnstate(self.b1))

    self.hl = QHBoxLayout()
    self.hl.addWidget(self.b1)
    self.hl.addWidget(self.sp)
    self.hl.addWidget(self.lbl)
    self.hl.addWidget(self.btn)
    self.setLayout(self.hl)

  def slider_moved(self):
    self.lbl.setText(str(self.sp.value()))

  def button_clicked(self, btn: QPushButton):
    if btn == self.btn:
      val = self.sp.value()*0.001
      self.cadmsh.edge_length = val
      self.cadmsh.remesh()
      if self.func_updated is not None:
        self.func_updated()

  def btnstate(self, b):
    if b == self.b1:
      self.cadmsh.is_sync_mesh = b.isChecked()


class QUI_ValueSlider(QWidget):
  def __init__(self,txt:str):
    super(QUI_ValueSlider, self).__init__()

    self.lbl1a = QLabel(txt+":")

    self.sp1 = QSlider(Qt.Horizontal)
    self.sp1.setMinimum(10)
    self.sp1.setMaximum(100)
    self.sp1.setValue(50)

    self.lbl1b = QLabel("50", self)

    self.hl1 = QHBoxLayout()
    self.hl1.addWidget(self.lbl1a)
    self.hl1.addWidget(self.sp1)
    self.hl1.addWidget(self.lbl1b)
    self.hl1.addStretch()
    self.setLayout(self.hl1)

class QUI_FEMPoisson(QWidget):
  def __init__(self,fem:delfem2.FEM_Poisson):
    super(QUI_FEMPoisson, self).__init__()
    self.fem = fem
    self.func_updated = None

#    self.b1 = QCheckBox("Sync Mesh To CAD")
#    self.b1.setChecked(True)
#    self.b1.stateChanged.connect(lambda: self.btnstate(self.b1))

    self.btn = QPushButton('solve', self)
    self.btn.clicked.connect(lambda: self.button_clicked(self.btn))

    self.hl0 = QHBoxLayout()
    self.hl0.addWidget(self.btn)

    ####

    self.vs1 = QUI_ValueSlider("alpha")
    self.vs1.sp1.sliderMoved.connect(lambda: self.slider_moved(self.vs1.sp1))

    ####

    self.vl = QVBoxLayout()
    self.vl.addLayout(self.hl0)
    self.vl.addWidget(self.vs1)
    self.setLayout(self.vl)

  def slider_moved(self,sp):
    if sp == self.vs1.sp1:
      self.vs1.lbl1b.setText(str(self.vs1.sp1.value()))
      self.fem.alpha = self.vs1.sp1.value() * 0.01

  def button_clicked(self, btn: QPushButton):
    if btn == self.btn:
      self.fem.solve()
      if self.func_updated is not None:
        self.func_updated()



class QUI_FEMSolidLinearStatic(QWidget):
  def __init__(self,fem:delfem2.FEM_SolidLinearStatic):
    super(QUI_FEMSolidLinearStatic, self).__init__()

    self.fem = fem
    self.func_updated = None

    self.btn = QPushButton('solve', self)
    self.btn.clicked.connect(lambda: self.button_clicked(self.btn))

    self.cb1 = QCheckBox("Sync FEM to Params")
    self.cb1.setChecked(True)

    self.hl0 = QHBoxLayout()
    self.hl0.addWidget(self.cb1)
    self.hl0.addStretch()
    self.hl0.addWidget(self.btn)

    ####

    self.vs1 = QUI_ValueSlider("myu")
    self.vs1.sp1.valueChanged.connect(lambda: self.slider_moved(self.vs1.sp1))

    self.vs2 = QUI_ValueSlider("lambda")
    self.vs2.sp1.valueChanged.connect(lambda: self.slider_moved(self.vs2.sp1))

    self.vs3 = QUI_ValueSlider("rho")
    self.vs3.sp1.valueChanged.connect(lambda: self.slider_moved(self.vs3.sp1))

    ####

    self.vl = QVBoxLayout()
    self.vl.addLayout(self.hl0)
    self.vl.addWidget(self.vs1)
    self.vl.addWidget(self.vs2)
    self.vl.addWidget(self.vs3)
    self.setLayout(self.vl)


  def slider_moved(self, sp: QSlider):
    if sp == self.vs1.sp1:
      self.vs1.lbl1b.setText(str(self.vs1.sp1.value()))
      self.fem.myu = self.vs1.sp1.value() * 0.01
    if sp == self.vs2.sp1:
      self.vs2.lbl1b.setText(str(self.vs2.sp1.value()))
      self.fem.lmd = self.vs2.sp1.value() * 0.01
    if sp == self.vs3.sp1:
      self.vs3.lbl1b.setText(str(self.vs3.sp1.value()))
      self.fem.rho = self.vs3.sp1.value() * 0.01
    if self.cb1.isChecked():
      self.fem.solve()
      if self.func_updated is not None:
        self.func_updated()

  def button_clicked(self, btn: QPushButton):
    if btn == self.btn:
      self.fem.solve()
      if self.func_updated is not None:
        self.func_updated()
