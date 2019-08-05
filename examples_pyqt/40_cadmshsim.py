####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import sys, math

from PySide2.QtCore import Qt, QSize, Signal
from PySide2.QtWidgets import QApplication, QWidget, QSizePolicy
from PySide2.QtWidgets import QVBoxLayout, QHBoxLayout, QCheckBox, QButtonGroup, QRadioButton

import PyDelFEM2 as dfm2
import PyDelFEM2.gl
import PyDelFEM2.qt
import PyDelFEM2.cadmshsimvis


class Window_Poisson(dfm2.qt.QW_CadMshFem):
  def __init__(self):
    super(Window_Poisson, self).__init__()
    self.setWindowTitle("CAD_Mesh_Poisson")

    self.cadmsh = dfm2.cadmshsimvis.CadMesh2D_FEMPoisson(edge_length=0.05)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])
    self.cadmsh.remesh()

    self.ui_fem = dfm2.qt.QW_FEMPoisson(self.cadmsh.fem)

    super().init_UI()


class Window_Diffuse(dfm2.qt.QW_CadMshFem):
  def __init__(self):
    super(Window_Diffuse, self).__init__()

    self.cadmsh = dfm2.cadmshsimvis.CadMesh2D_FEMDiffuse(edge_length=0.05)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])
    self.cadmsh.remesh()

    self.ui_fem = dfm2.qt.QW_FEMDiffuse(self.cadmsh.fem)

    self.setWindowTitle("CAD_Mesh_Diffuse")

    super().init_UI()


class Window_SolidLinearStatic(dfm2.qt.QW_CadMshFem):
  def __init__(self):
    super(Window_SolidLinearStatic, self).__init__()

    self.setWindowTitle("CAD_Mesh_SolidLinearStatic")

    self.cadmsh = dfm2.cadmshsimvis.CadMesh2D_FEMSolidLinearStatic(edge_length=0.05)
    self.cadmsh.fem.param_gravity_y = -0.1
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])
    self.cadmsh.remesh()

    self.ui_fem = dfm2.qt.QW_FEMSolidLinearStatic(self.cadmsh.fem)

    super().init_UI()


class Window_SolidLinearEigen(dfm2.qt.QW_CadMshFem):
  def __init__(self):
    super(Window_SolidLinearEigen, self).__init__()
    self.setWindowTitle("CAD_Mesh_SolidLinearEigen")

    self.cadmsh = dfm2.cadmshsimvis.CadMesh2D_FEMSolidLinearEigen(edge_length=0.05)
    self.cadmsh.add_polygon([-1, -0.2, +1, -0.2, +1, +0.2, -1, +0.2])
    self.cadmsh.remesh()

    self.ui_fem = dfm2.qt.QW_FEMSolidLinearEigen(self.cadmsh.fem)

    super().init_UI()


class Window_Pbd2D(dfm2.qt.QW_CadMshFem):
  def __init__(self):
    super(Window_Pbd2D, self).__init__()

    self.setWindowTitle("CAD_Mesh_Pbd2d")

    self.cadmsh = dfm2.cadmshsimvis.CadMesh2D_PBD(edge_length=0.1)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, -1, +1])
    self.cadmsh.remesh()

    self.ui_fem = dfm2.qt.QW_PBD(self.cadmsh.pbd)

    super().init_UI()


class Window_PbdCloth(dfm2.qt.QW_CadMshFem):
  def __init__(self):
    super(Window_PbdCloth, self).__init__()

    self.setWindowTitle("CAD_Mesh_Pbd2d")

    self.cadmsh = dfm2.cadmshsimvis.CadMesh2D_PBDCloth(edge_length=0.1)
    self.cadmsh.add_polygon([-1, -1, +1, -1, +1, +1, +0.8,+1, -0.8,+1, -1, +1])
    self.cadmsh.listIndEdge_Fix = [2,4]
    self.cadmsh.remesh()

    self.cadmsh.pbd.param_gravity_y = -0.1
    self.cadmsh.pbd.param_gravity_z = -0.001
    self.cadmsh.pbd.dt = 0.1

    self.ui_fem = dfm2.qt.QW_PBDCloth(self.cadmsh.pbd)

    super().init_UI()


if __name__ == '__main__':
  app = QApplication(sys.argv)

  window = Window_Poisson()
  window.show()
  app.exec_()

  window = Window_Diffuse()
  window.show()
  app.exec_()

  window = Window_SolidLinearStatic()
  window.show()
  app.exec_()

  window = Window_SolidLinearEigen()
  window.show()
  app.exec_()

  window = Window_Pbd2D()
  window.show()
  app.exec_()

  window = Window_PbdCloth()
  window.show()
  app.exec_()