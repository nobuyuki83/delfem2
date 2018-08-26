from OpenGL.GL import *
from OpenGL.GLUT import *

import sys
sys.path.append("../python")
import dfm2

def draw_func():
  glEnable(GL_LIGHTING)
  glutSolidTeapot(0.5)

def main():
  win = dfm2.WindowGLFW(1.0);
  dfm2.setSomeLighting()
  win.draw_loop(draw_func)

if __name__ == "__main__":  
  main()
