import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

import numpy

gltf = PyDelFEM2.CppGLTF()
#gltf.read("../test_inputs/RiggedFigure.glb")
gltf.read("../test_inputs/CesiumMan.glb")
#gltf.print()

np_pos0,np_elm,np_rigw,np_rigj = dfm2.CppGLTF_GetMeshInfo(gltf,0,0)
np_pos = np_pos0.copy()
msh = dfm2.Mesh(np_pos,np_elm,dfm2.TRI)
dfm2.gl.glfw.winDraw3d([msh],winsize=(400,300))

bone_array = dfm2.CppGLTF_GetBones(gltf,0)
bone_array.set_rotation_bryant(0, [-3.1415*0.5, 0.0, 0.0] )
bone_array.set_translation(0, [0.0, 0.0, +0.2])
dfm2.update_rig_skin(np_pos,
                    np_pos0,np_elm,bone_array,np_rigw,np_rigj)
msh = dfm2.Mesh(np_pos,np_elm,dfm2.TRI)
axis = dfm2.gl.AxisXYZ(1)
dfm2.gl.glfw.winDraw3d([msh,axis],winsize=(400,300))