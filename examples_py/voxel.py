from OpenGL.GL import *

import sys
sys.path.append("../python")
import dfm2

mshelm = None

def display():
  glEnable(GL_LIGHTING)
  mshelm.drawFace_elemWiseNorm()
  glDisable(GL_LIGHTING)
  mshelm.drawEdge()

  glDisable(GL_LIGHTING)
  glLineWidth(3)  
  dfm2.draw_axis(size=0.5)

def main():
  global mshelm

  voxelgrid = dfm2.VoxelGrid()
  voxelgrid.add(0, 0, 0)
  voxelgrid.add(1, 0, 0)
  voxelgrid.add(0, 1, 0)

  mshelm = dfm2.meshQuad3d_voxelGrid(voxelgrid)
  mshelm = mshelm.subdiv()
  mshelm = mshelm.subdiv()  
  mshelm = mshelm.subdiv()  

  win = dfm2.WindowGLFW(1.0)
  dfm2.setSomeLighting()
  glEnable(GL_POLYGON_OFFSET_FILL );
  glPolygonOffset( 1.1, 4.0 );
  win.draw_loop(display)

if __name__ == "__main__":
  main()