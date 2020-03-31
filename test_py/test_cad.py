####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import pytest
import numpy, random, os, sys, math
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

os.chdir(os.path.dirname(os.path.abspath(__file__))) # for python3 setup.py test


class Test_CppCad2D():
  def test1(self,request):
    is_gl = request.config.getoption('--is_gl') == "true"

    ccad = dfm2.CppCad2D()
    ccad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
    assert ccad.check() 
    if is_gl:
        dfm2.gl.glfw.winDraw3d([ccad],nframe=20)

    ccad.add_vtx_face(0.0, 0.0, 0)
    assert ccad.check()
    if is_gl:
        dfm2.gl.glfw.winDraw3d([ccad], nframe=20)

    ccad.add_vtx_edge(0.0,0.8, 2)
    assert ccad.check()
    if is_gl:
        dfm2.gl.glfw.winDraw3d([ccad],nframe=20)

    ccad.set_edge_type(0, dfm2.CAD_EDGE_GEOM_BEZIER_CUBIC, [0.2, 0.3, -0.2, 0.3])
    assert ccad.check()
    if is_gl:
        dfm2.gl.glfw.winDraw3d([ccad], nframe=20)

    ccad.set_edge_type(0, dfm2.CAD_EDGE_GEOM_LINE, [])
    assert ccad.check()
    if is_gl:
        dfm2.gl.glfw.winDraw3d([ccad], nframe=20)


class Test_Cad2D():
  def test1(self):
    cad = dfm2.Cad2D()
    cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
    assert cad.ccad.ind_vtx_face(0) == [0,1,2,3]
    assert cad.ccad.ind_edge_face(0) == [(0,True),(1,True),(2,True),(3,True)]
    assert cad.ccad.check()
    ####
    cad.add_vtx_edge(0,[0.0,0.0])
    assert cad.ccad.check()
    assert cad.ccad.ind_vtx_face(0) == [0,4,1,2,3]
    assert cad.ccad.ind_edge_face(0) == [(0,True),(4,True),(1,True),(2,True),(3,True)]

  def test2(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
    assert cad.ccad.check()
    ####
    cad.pick(-1,-1, 2.0)
    assert cad.ivtx_picked() == 0 
    assert cad.iedge_picked() == -1 
    assert cad.iface_picked() == -1 
    ####
    cad.pick(+0,-1, 2.0)
    assert cad.ivtx_picked() == -1 
    assert cad.iedge_picked() == 0 
    assert cad.iface_picked() == -1 
    ####
    cad.pick(-0.1432143,-0.1653542343, 2.0)
    assert cad.ivtx_picked() == -1 
    assert cad.iedge_picked() == -1 
    assert cad.iface_picked() == 0 

  def test3(self):
    cad = dfm2.Cad2D()
    cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
    mesher = dfm2.Mesher_Cad2D(edge_length=0.05)
    dmsh = mesher.meshing(cad)
    assert dmsh.np_pos.shape[1] == 2
    np_xy_bound = numpy.array(cad.ccad.xy_vtxctrl_face(0)).reshape([-1, 2])
    W = dfm2.mvc(dmsh.np_pos, np_xy_bound)
    assert W.ndim == 2
    assert W.shape[0] == dmsh.np_pos.shape[0]
    assert W.shape[1] == len(cad.ccad.xy_vtxctrl_face(0))/2
    assert numpy.linalg.norm(W.sum(axis=1)-numpy.ones((W.shape[0])))<1.0e-3


  def test_svg(self,request):
    is_gl = request.config.getoption('--is_gl') == "true"
    cad = dfm2.Cad2D()
    for itr in range(3):
      cad.clear()
      cad.import_svg("../test_inputs/shape" + str(itr) + ".svg")
      if is_gl:
        dfm2.gl.glfw.winDraw3d([cad],nframe=20)
