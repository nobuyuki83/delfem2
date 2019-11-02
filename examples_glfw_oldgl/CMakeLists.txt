cmake_minimum_required(VERSION 2.8)

project(examples_glfw_oldgl)

# simple examples without delfem2 
add_subdirectory(000_openwin)

# basic funcs
add_subdirectory(100_nav3d)
add_subdirectory(101_subdiv)
add_subdirectory(102_convhull)
add_subdirectory(103_noise2d)
add_subdirectory(104_noise3d)
add_subdirectory(105_meshslice)
add_subdirectory(106_sphericalharmonics)
add_subdirectory(107_paramgeo2d)
add_subdirectory(108_paramgeo3d)

# meshing
add_subdirectory(201_dyntri2d_triangulation)
add_subdirectory(202_dyntri2d_remesh)
add_subdirectory(203_dyntri3d_edgecollapse)
add_subdirectory(204_dyntet_tetrahedralization)