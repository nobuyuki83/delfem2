####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys, math, numpy
import OpenGL.GL as gl

from typing import List

from PySide2.QtCore import QPoint, QSize, Qt, QEvent, QTimer, QObject, Signal
from PySide2.QtGui import QColor
from PySide2.QtWidgets import QOpenGLWidget, QMenu, QWidget, QSizePolicy
from PySide2.QtWidgets import \
  QPushButton, QLabel, QSlider, QHBoxLayout, QCheckBox, QVBoxLayout, QGridLayout, \
  QLabel, QRadioButton,\
  QButtonGroup

import PyDelFEM2 as dfm2
import PyDelFEM2.gl

def setClearColor(c:QColor):
  gl.glClearColor(c.redF(), c.greenF(), c.blueF(), c.alphaF())

def setColor(c:QColor):
  gl.glColor4f(c.redF(), c.greenF(), c.blueF(), c.alphaF())


class NavigationPyQt(QObject):
  updated = Signal()
  def __init__(self,view_height):
    super(NavigationPyQt,self).__init__()
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
    self.updated.emit()

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
    self.updated.emit()

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


class QGLW_Nav3D(QOpenGLWidget):
  def __init__(self, parent=None, nav=None):
    super(QGLW_Nav3D, self).__init__(parent)
    if nav is None:
      self.nav = NavigationPyQt(view_height=1.0)
    else:
      self.nav = nav

  def minimumSizeHint(self):
    return QSize(300, 300)

  def sizeHint(self):
    return QSize(600, 600)

  def initializeGL(self):
    gl.glClearColor(0.8, 0.8, 1.0, 1.0)
    gl.glShadeModel(gl.GL_FLAT)
    gl.glEnable(gl.GL_DEPTH_TEST)
    gl.glEnable(gl.GL_CULL_FACE)
    gl.glDisable(gl.GL_LIGHTING)
    gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
    gl.glPolygonOffset(1.1, 4.0 )
    dfm2.gl.setSomeLighting()

  def paintGL(self):
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    self.nav.camera.set_gl_camera()

  def resizeGL(self, width, height):
    gl.glViewport(0,0,width,height)

  def mousePressEvent(self, event):
    self.nav.mouse(event, self.frameSize())

  def mouseMoveEvent(self, event):
    self.nav.motion(event, self.frameSize())
    self.update()

  def wheelEvent(self,event):
    v = 0.1*(event.angleDelta().y()/120.0)
    self.nav.camera.scale *= math.exp(v)
    self.nav.updated.emit()
    self.update()




class QGLW_Mesh(QGLW_Nav3D):
  def __init__(self):
    super(QGLW_Mesh, self).__init__()
    self.msh = None

  def initializeGL(self):
    super().initializeGL()
    self.buff_face = dfm2.GLBufferMesh(self.msh, is_normal=True)

  def paintGL(self):
    super().paintGL()
    if self.buff_face is None:
      self.msh.draw()
    else:
      gl.glEnable(gl.GL_LIGHTING)
      self.buff_face.draw()



class QGLW_Cad2D(QOpenGLWidget):
  updated_cadmshfem = Signal()

  def __init__(self, parent=None):
    super(QGLW_Cad2D, self).__init__(parent)
    self.cadobj = None
    self.nav = dfm2.qt.NavigationPyQt(view_height=2.0)
    self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

  def minimumSizeHint(self):
    return QSize(300, 300)

  def sizeHint(self):
    return QSize(2000, 2000)

  def initializeGL(self):
    print(dfm2.gl.getOpenglInfo())
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
    self.updated_cadmshfem.emit()

  def mouseReleaseEvent(self, event):
#    self.cadobj.clean_picked()
#    self.update()
    pass

  def mouseMoveEvent(self, event):
    self.makeCurrent()
    self.nav.motion(event, self.frameSize())
    src0,src1,dir = self.nav.motion_src_dir()
    self.cadobj.motion(src0, src1, dir)
    self.update()
    self.updated_cadmshfem.emit()

  def contextMenuEvent(self, event):
    src,dir = self.nav.mouse_src_dir()
    ###
    menu = QMenu(self)
    actionDelVtx = None
    actionAddVtx = None
    actionEdgeBezier = None
    actionAddSquare = None
    if self.cadobj.ivtx_picked() != -1:
      actionDelVtx = menu.addAction("delete vtx")
    elif self.cadobj.iedge_picked() != -1:
      actionAddVtx = menu.addAction("add vtx")
      edge_type = self.cadobj.edge_type( self.cadobj.iedge_picked() )
      if edge_type ==0:
        actionEdgeBezier = menu.addAction("set bezier")
    else:
      actionAddSquare = menu.addAction("add square")
    ####
    action = menu.exec_(self.mapToGlobal(event.pos()))
    if action == actionAddVtx != None:
      self.cadobj.add_vtx_edge(iedge=self.cadobj.ccad.iedge_picked,
                               pos=[src[0], src[1]] )
    if action == actionAddSquare != None:
      x0,y0 = src[0],src[1]
      self.cadobj.add_polygon([x0-1, y0-1, x0+1, y0-1, x0+1, y0+1, x0-1, y0+1])
    if action == actionEdgeBezier != None:
      self.cadobj.set_edge_type( self.cadobj.iedge_picked(), 1, [0.2,0.3,-0.2,0.3] )
    self.update()
    self.updated_cadmshfem.emit()

  def wheelEvent(self,event):
    v = 0.1*(event.angleDelta().y()/120.0)
    self.nav.camera.scale *= math.exp(v)
    self.update()





class QW_NumWin(QWidget):
  changed = Signal(int)

  def __init__(self, numwin:int):
    super(QW_NumWin,self).__init__()
    r1=QRadioButton("one window")
    r1.toggled.connect(self.checked)
    r2=QRadioButton("two window")
    r2.toggled.connect(self.checked)
    button_group_numwin = QButtonGroup()
    button_group_numwin.addButton(r1)
    button_group_numwin.addButton(r2)

    if numwin == 1:
      r1.setChecked(True)
    else:
      r2.setChecked(True)

    self.layout_numwin = QHBoxLayout()
    self.layout_numwin.addWidget(r1, alignment=Qt.AlignLeft)
    self.layout_numwin.addWidget(r2, alignment=Qt.AlignLeft)
    self.setLayout(self.layout_numwin)

  def checked(self):
    ch = self.sender()
    if ch.text() == "one window" and ch.isChecked():
      self.changed.emit(1)
    elif ch.text() == "two window" and ch.isChecked():
      self.changed.emit(2)


class QW_ValudSlider(QWidget):
  def __init__(self,name:str):
    super(QW_ValudSlider, self).__init__()

    self.name = name

    self.lbl1a = QLabel(name+":")

    self.sp1 = QSlider(Qt.Horizontal)
    self.sp1.setMinimum(10)
    self.sp1.setMaximum(100)
    self.sp1.setValue(50)

    self.lbl1b = QLabel("50", self)

    self.hl1 = QHBoxLayout()
    self.hl1.addWidget(self.lbl1a, alignment=Qt.AlignLeft)
    self.hl1.addWidget(self.sp1,   alignment=Qt.AlignLeft)
    self.hl1.addWidget(self.lbl1b, alignment=Qt.AlignLeft)
    self.hl1.addStretch()
    self.setLayout(self.hl1)


class QW_FemParams(QWidget):
  valueChanged = Signal()

  def __init__(self,list_name:List[str],fem):
    super(QW_FemParams, self).__init__()

    self.list_name = list_name
    self.fem = fem
    self.grid = QGridLayout()

    self.list_sldr = []
    for iname,name in enumerate(list_name):
      lbl1a = QLabel(name+":")
      ##
      sp1 = QSlider(Qt.Horizontal)
      sp1.setMinimum(10)
      sp1.setMaximum(100)
      sp1.setValue(50)
      sp1.valueChanged.connect(self.sliderValueChanged)
      self.list_sldr.append(sp1)
      ##
      lbl1b = QLabel("50", self)
      ##
      self.grid.addWidget(lbl1a, iname,0)
      self.grid.addWidget(sp1,   iname,1)
      self.grid.addWidget(lbl1b, iname,2)
    self.setLayout(self.grid)

  def sliderValueChanged(self,ival:int):
    slider = self.sender()
    ilist0 = self.list_sldr.index(slider)
    name = self.list_name[ilist0]
    setattr(self.fem, "param_"+name, ival*0.01)
    self.valueChanged.emit()


class QW_MeshRes(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,cadmsh:dfm2.CadMesh2D):
    super(QW_MeshRes, self).__init__()
    self.cadmsh = cadmsh

    self.btn = QPushButton('remesh', self)
    self.btn.clicked.connect(lambda: self.button_clicked(self.btn))

    self.b1 = QCheckBox("Sync Mesh To CAD")
    self.b1.setChecked(True)
    self.b1.stateChanged.connect(lambda: self.btnstate(self.b1))

    self.vs = QW_ValudSlider("mesh density")
    self.vs.sp1.valueChanged.connect(self.slider_moved)

    self.hl = QHBoxLayout()
    self.hl.addWidget(self.b1, alignment=Qt.AlignLeft)
    self.hl.addWidget(self.btn,alignment=Qt.AlignLeft)

    self.vl = QVBoxLayout()
    self.vl.addLayout(self.hl)
    self.vl.addWidget(self.vs, alignment=Qt.AlignLeft)
    self.setLayout(self.vl)


  def slider_moved(self):
    self.vs.lbl1b.setText(str(self.vs.sp1.value()))

  def button_clicked(self, btn: QPushButton):
    if btn == self.btn:
      val = self.vs.sp1.value()*0.001
      self.cadmsh.edge_length = val
      self.cadmsh.remesh()
      self.updated_cadmshfem.emit()

  def btnstate(self, b):
    if b == self.b1:
      self.cadmsh.is_sync_mesh = b.isChecked()


class QW_SolveParam(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,fem):
    super(QW_SolveParam,self).__init__()

    self.fem = fem

    self.btn = QPushButton('solve', self)
    self.btn.pressed.connect(self.button_clicked)

    self.cb1 = QCheckBox("Sync FEM to Params")
    self.cb1.setChecked(True)

    self.hl0 = QHBoxLayout()
    self.hl0.addWidget(self.cb1, alignment=Qt.AlignLeft)
    self.hl0.addStretch()
    self.hl0.addWidget(self.btn, alignment=Qt.AlignLeft)
    self.setLayout(self.hl0)

  def button_clicked(self):
    self.fem.solve()
    self.updated_cadmshfem.emit()


class QW_AnimCntrl(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,fem):
    super(QW_AnimCntrl,self).__init__()

    self.fem = fem

    self.tmr_stepTime = QTimer(self)
    self.tmr_stepTime.setSingleShot(False)
    self.tmr_stepTime.timeout.connect(self.stepTime)

    self.btn_Initialize = QPushButton("initialize")

    self.btn_Animate = QPushButton("animate")
    self.btn_Animate.setCheckable(True)
    self.btn_Animate.toggled.connect(self.btn_Animate_toggled)

    self.btn_Initialize.pressed.connect(self.btn_Initialize_pressed)

    self.hl = QHBoxLayout()
    self.hl.addWidget(self.btn_Initialize)
    self.hl.addWidget(self.btn_Animate)
    self.setLayout(self.hl)

  def btn_Animate_toggled(self):
    if self.btn_Animate.isChecked():
      self.tmr_stepTime.start(30)
    else:
      self.tmr_stepTime.stop()

  def btn_Initialize_pressed(self):
    self.fem.initialize()
    self.updated_cadmshfem.emit()

  def stepTime(self):
    self.fem.step_time()
    self.updated_cadmshfem.emit()

##########################################################################################


class QW_FEMPoisson(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,fem:dfm2.FEM_Poisson):
    super(QW_FEMPoisson, self).__init__()
    self.fem = fem
    
    self.cb1 = QCheckBox("Sync FEM to Params")
    self.cb1.setChecked(True)
    self.cb1.stateChanged.connect(lambda: self.btnstate(self.b1))

    self.qui_sp = QW_SolveParam(self.fem)
    self.qui_sp.updated_cadmshfem.connect(lambda: self.updated_cadmshfem.emit())

    ####

    self.vs1 = QW_FemParams(["alpha","source"],self.fem)
    self.vs1.valueChanged.connect(self.fem_param_changed)

    ####

    self.vl = QVBoxLayout()
    self.vl.addWidget(self.qui_sp, alignment=Qt.AlignLeft)
    self.vl.addWidget(self.vs1,    alignment=Qt.AlignLeft)
    self.vl.addStretch()
    self.setLayout(self.vl)


  def fem_param_changed(self):
    if self.qui_sp.cb1.isChecked():
      self.fem.solve()
      self.updated_cadmshfem.emit()




class QW_FEMSolidLinearStatic(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,fem:dfm2.FEM_SolidLinearStatic):
    super(QW_FEMSolidLinearStatic, self).__init__()

    self.fem = fem

    ####
    self.qui_sp = QW_SolveParam(self.fem)
    self.qui_sp.updated_cadmshfem.connect(lambda: self.updated_cadmshfem.emit())
    ####
    self.vs1 = QW_FemParams(["myu","lambda","rho"],self.fem)
    self.vs1.valueChanged.connect(self.fem_param_changed)
    ####

    self.vl = QVBoxLayout()
    self.vl.addWidget(self.qui_sp,  alignment=Qt.AlignLeft)
    self.vl.addWidget(self.vs1,alignment=Qt.AlignLeft)
    self.setLayout(self.vl)

  def fem_param_changed(self):
    if self.qui_sp.cb1.isChecked():
      self.fem.solve()
      self.updated_cadmshfem.emit()


class QW_FEMSolidLinearEigen(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,fem:dfm2.FEM_SolidLinearEigen):
    super(QW_FEMSolidLinearEigen, self).__init__()
    ####
    self.fem = fem
    ####
    self.qui_animctrl = QW_AnimCntrl(self.fem)
    self.qui_animctrl.updated_cadmshfem.connect(lambda: self.updated_cadmshfem.emit())
    ####
    self.vs1 = QW_FemParams(["myu","lambda","rho"],self.fem)
    self.vs1.valueChanged.connect(lambda: self.updated_cadmshfem.emit())
    ####
    self.vl = QVBoxLayout()
    self.vl.addWidget(self.qui_animctrl,  alignment=Qt.AlignLeft)
    self.vl.addWidget(self.vs1,alignment=Qt.AlignLeft)
    self.setLayout(self.vl)


class QW_FEMDiffuse(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,fem:dfm2.FEM_Diffuse):
    super(QW_FEMDiffuse, self).__init__()
    ####
    self.fem = fem
    ####
    self.qui_animctrl = QW_AnimCntrl(self.fem)
    self.qui_animctrl.updated_cadmshfem.connect(lambda: self.updated_cadmshfem.emit())
    ####
    self.vs1 = QW_FemParams(["alpha","source"],self.fem)
    self.vs1.valueChanged.connect(lambda: self.updated_cadmshfem.emit())
    ####
    self.vl = QVBoxLayout()
    self.vl.addWidget(self.qui_animctrl, alignment=Qt.AlignLeft)
    self.vl.addWidget(self.vs1,alignment=Qt.AlignLeft)
    self.setLayout(self.vl)



#####################################################


class QW_PBD(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,fem:dfm2.PBD):
    super(QW_PBD, self).__init__()

    self.fem = fem

    self.qui_sp = QW_SolveParam(self.fem)
    self.qui_sp.btn.clicked.connect(lambda: self.button_clicked(self.qui_sp.btn))
    self.qui_sp.updated_cadmshfem.connect(lambda: self.updated_cadmshfem.emit())

    self.qui_InitAnim = QW_AnimCntrl(self.fem)
    self.qui_InitAnim.updated_cadmshfem.connect(lambda: self.updated_cadmshfem.emit())


    self.vl = QVBoxLayout()
    self.vl.addWidget(self.qui_InitAnim)
#    self.vl.addWidget(self.vs1,alignment=Qt.AlignLeft)
    self.setLayout(self.vl)

'''    
  def fem_param_changed(self):
    if self.qui_sp.cb1.isChecked():
      self.fem.solve()
      self.updated_cadmshfem.emit()

  def button_clicked(self, btn: QPushButton):
    if btn == self.qui_sp.btn:
      self.fem.solve()
      self.updated_cadmshfem.emit()
'''



class QW_PBDCloth(QWidget):
  updated_cadmshfem = Signal()

  def __init__(self,fem:dfm2.PBD_Cloth):
    super(QW_PBDCloth, self).__init__()

    self.fem = fem

    ####

    self.qui_animctrl = QW_AnimCntrl(self.fem)
    self.qui_animctrl.updated_cadmshfem.connect(lambda: self.updated_cadmshfem.emit())

#    self.vs1 = QW_FemParams(["alpha","source"],self.fem)
#    self.vs1.valueChanged.connect(self.fem_param_changed)

    ####

    self.vl = QVBoxLayout()
    self.vl.addWidget(self.qui_animctrl, alignment=Qt.AlignLeft)
#    self.vl.addWidget(self.vs1,alignment=Qt.AlignLeft)
    self.setLayout(self.vl)

  def fem_param_changed(self):
    self.updated_cadmshfem.emit()


#####################################################


class QW_CadMshFem(QWidget):
  def __init__(self):
    super(QW_CadMshFem, self).__init__()
    self.setWindowTitle("CadMshFem")

  def sizeHint(self):
    return QSize(1200, 600)

  def init_UI(self):
    self.glWidget0 = dfm2.qt.QGLW_Cad2D()
    self.glWidget0.cadobj = self.cadmsh
    self.glWidget0.updated_cadmshfem.connect(self.updated_cadmshfem)

    self.glWidget1 = dfm2.qt.QGLW_Cad2D()
    self.glWidget1.cadobj = self.cadmsh
    self.glWidget1.updated_cadmshfem.connect(self.updated_cadmshfem)

    self.ui_numwin = dfm2.qt.QW_NumWin(1)
    self.ui_numwin.changed.connect(self.numwin_changed)
    self.numwin_changed(1)

    self.ui_meshres = dfm2.qt.QW_MeshRes(self.cadmsh)
    self.ui_meshres.updated_cadmshfem.connect(self.updated_cadmshfem)

    self.layout_param = QVBoxLayout()
    self.layout_param.addWidget(self.ui_numwin,  alignment=Qt.AlignLeft)
    self.layout_param.addWidget(self.ui_meshres, alignment=Qt.AlignLeft)
    if hasattr(self,"ui_fem"):
      self.layout_param.addWidget(self.ui_fem,     alignment=Qt.AlignLeft)
      self.ui_fem.updated_cadmshfem.connect(self.updated_cadmshfem)
    self.layout_param.addStretch()

    self.mainLayout = QHBoxLayout()
    self.mainLayout.addWidget(self.glWidget0, alignment=Qt.AlignLeft, stretch=2)
    self.mainLayout.addWidget(self.glWidget1, alignment=Qt.AlignLeft, stretch=2)
    self.mainLayout.addLayout(self.layout_param, stretch=0)
    self.setLayout(self.mainLayout)

  def updated_cadmshfem(self):
    self.glWidget0.update()
    self.glWidget1.update()

  def keyPressEvent(self, event):
    if event.text() == 'q':
      self.close()

  def numwin_changed(self,numwin:int):
    if numwin == 1:
      self.glWidget1.hide()
    elif numwin == 2:
      self.glWidget1.show()

