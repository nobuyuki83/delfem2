from OpenGL.GL import *
from OpenGL.GLUT import *

import sys
sys.path.append("../python")
import dfm2

win = None

def draw_func():
  glEnable(GL_LIGHTING)
  glutSolidTeapot(0.5)

def main():
  draw = dfm2.glut.GlutWindow(1.0);
  dfm2.dfm2.setSomeLighting()
  draw.draw_loop(draw_func)

if __name__ == "__main__":  
  main()
