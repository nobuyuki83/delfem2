import sys
sys.path.append("../module_py")
import dfm2

def main():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])

  cadmsh = dfm2.CadMeshLinked(cad,elen=0.1)

  dfm2.winDraw3d([cadmsh])

if __name__ == "__main__":
  main()