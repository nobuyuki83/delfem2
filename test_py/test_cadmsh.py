####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

#import unittest,
import pytest
import numpy, random, os, sys, math
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

os.chdir(os.path.dirname(os.path.abspath(__file__))) # for python3 setup.py test


class Test_MathExpression():
  def test1(self):
    mee = dfm2.MathExpressionEvaluator()
    mee.set_key("x",3)
    mee.set_expression("x+3")
    mee.set_key("x",6)
    assert math.fabs(mee.eval()-9.0)<1.0e-30


class Test_Mesh():

  def test1(self,request):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_2k.ply")
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)


  def test2(self,request):
    msh = dfm2.Mesh()
    msh.read("../test_inputs/bunny_1k.obj")
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test3(self,request):
    msh = dfm2.Mesh()
    msh.set_grid((32,16))
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

  def test4(self,request):
    voxelgrid = dfm2.VoxelGrid()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = voxelgrid.mesh_quad()
    msh = msh.subdiv()
    msh = msh.subdiv()
    msh = msh.subdiv()
    assert msh is not None
    if request.config.getoption('--is_gl') == "true":    
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)

'''
  def test5(self):
    voxelgrid = dfm2.VoxelGrid()
    voxelgrid.add(0, 0, 0)
    voxelgrid.add(1, 0, 0)
    voxelgrid.add(0, 1, 0)
    msh = voxelgrid.mesh_hex()
    msh = msh.subdiv()
    self.assertIsNot(msh,None)
    if self.is_gl:
      dfm2.gl.glfw.winDraw3d([msh],nframe=20)
'''      
