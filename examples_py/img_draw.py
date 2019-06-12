import matplotlib.pyplot as plt
import sys
sys.path.append("../module_py")
import delfem2 as dfm2
import delfem2.glfw

msh = dfm2.Mesh()
msh.read("../test_inputs/bunny_2k.ply")
img = dfm2.glfw.imgDraw3d([msh],winsize=(400,300))
plt.imshow(img)
plt.show()
