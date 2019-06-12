import sys
sys.path.append("../module_py")
import dfm2
import dfm2.glfw

loop = [[0,0, 1,0, 1,1, 0,1]]
out = dfm2.triangulation(loop,edge_length=0.05)
msh = dfm2.Mesh()
msh.meshtri2d(out[0],out[1])
axis = dfm2.AxisXYZ(1)
dfm2.glfw.winDraw3d([msh, axis])

