import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw


gltf = PyDelFEM2.CppGLTF()
gltf.read("../test_inputs/RiggedFigure.glb")
gltf.print()
