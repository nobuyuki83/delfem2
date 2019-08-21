import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw


gltf = PyDelFEM2.CppGLTF()
gltf.read("../test_inputs/RiggedFigure.glb")
#gltf.print()
np_pos,np_elm,np_rigw,np_rigj = dfm2.CppGLTF_GetMeshInfo(gltf,0,0)
msh = dfm2.Mesh(np_pos,np_elm,dfm2.TRI)
dfm2.gl.glfw.winDraw3d([msh],winsize=(400,300))
