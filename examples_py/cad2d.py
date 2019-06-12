import sys
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw

def main():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  dfm2.glfw.winDraw3d([cad])

if __name__ == "__main__":
  main()
