import sys
sys.path.append("../module_py")
import dfm2

def main():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])

  mesh = cad.mesh(0.1)
  dfm2.winDraw3d([cad,mesh])

if __name__ == "__main__":
  main()
