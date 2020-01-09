cmake_minimum_required(VERSION 2.8)

project(examples_glfw)

add_subdirectory(00_openwin)
add_subdirectory(01_drawrect)
add_subdirectory(02_nav3d)
add_subdirectory(03_texture)
add_subdirectory(04_render2texture)
add_subdirectory(05_OffScreenProjection)

add_subdirectory(10_imgui_sampledemo)
add_subdirectory(11_imgui_basic)

add_subdirectory(20_cad2d)

add_subdirectory(30_cadmsh2d)

add_subdirectory(40_femcloth)
add_subdirectory(41_fem2d_poisson)
add_subdirectory(42_fem2d_linearsolid)
