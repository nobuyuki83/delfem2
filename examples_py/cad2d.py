import sys
sys.path.append("../module_py")
import dfm2
import dfm2.glfw

def main():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
  dfm2.glfw.winDraw3d([cad])

if __name__ == "__main__":
  main()
