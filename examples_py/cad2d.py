import sys
sys.path.append("../module_py")
import dfm2

def main():
  cad = dfm2.Cad2D()
  cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])

#  xy,tri = dfm2.getMesh_cad(cad,0.1)
#  mesh = dfm2.Mesh(xy,tri)

  dfm2.winDraw3d([cad])


if __name__ == "__main__":
  main()
