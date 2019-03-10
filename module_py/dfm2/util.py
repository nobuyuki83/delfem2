import numpy
from .dfm2 import *



def normalize_rigmsh(rigmsh):
  aabb = rigmsh.aabb()
  center = aabb.center()
  rigmsh.translate(-center[0],-center[1],-center[2])
  rigmsh.scale(1.0/aabb.max_length())


class Field():
  def __init__(self,
               mesh: MeshElem,
               val: numpy.ndarray):
    self.mesh = mesh
    self.val = val
    self.draw_val_min = val.min()
    self.draw_val_max = val.max()
    self.color_map = ColorMap(self.draw_val_min,self.draw_val_max)

  def draw(self):
    draw_field(self.mesh,self.val,self.color_map)

