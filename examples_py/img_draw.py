import matplotlib.pyplot as plt
import sys
sys.path.append("../module_py")
import dfm2

msh = dfm2.Mesh("../test_inputs/bunny_2k.ply")
img = dfm2.imgDraw3d([msh],winsize=(400,300))
plt.imshow(img)
plt.show()
