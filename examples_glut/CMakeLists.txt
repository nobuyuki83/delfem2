cmake_minimum_required(VERSION 2.8)

project(examples_glut)

# opengl
add_subdirectory(00_viewer3d)
add_subdirectory(01_shader)

add_subdirectory(offscreen)
add_subdirectory(glbuffer)

# geometry
add_subdirectory(exponentialmap)
add_subdirectory(edit_cubes)
add_subdirectory(transform_handler)


# fem
add_subdirectory(rigidbody3d_rotation)

add_subdirectory(fem2d)
add_subdirectory(fem3d)
add_subdirectory(fem25d_eigen)
add_subdirectory(fem25d_stiffwarp)
add_subdirectory(femcloth_internal)
add_subdirectory(femcloth_selfcollision)
add_subdirectory(femcloth_dtri)
add_subdirectory(fempb_mitc3eigen)

# pbd
add_subdirectory(pbd_grid2d)
add_subdirectory(pbdcloth)
add_subdirectory(pbdcloth_cad)
add_subdirectory(pbdcloth_cad_gltf)
