from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

#msh_surf = dfm2.mesh_read("../test_inputs/bunny_2k.ply");
#msh_surf.scale_xyz(0.03)
#msh_tet = dfm2.mesh_isosurfaceStuffing(msh_surf.np_pos, msh_surf.np_elm)
#dfm2.winDraw3d([msh_surf],winsize=(400,300))

sdf = dfm2.SDF()
sdf.add_sphere(xyz=[0,0,0],rad=1)
np_xyz,np_tet = dfm2.isosurface(sdf.list_sdf)
print(np_xyz.shape,np_tet.shape)
msh = dfm2.Mesh(np_xyz,np_tet)

dfm2.winDraw3d([sdf,msh],winsize=(400,300))