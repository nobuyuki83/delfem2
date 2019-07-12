cmake_minimum_required(VERSION 2.8)

project(delfem2)

add_subdirectory(viewer3d_glut)
add_subdirectory(viewer3d_glfw)

add_subdirectory(noise2d)
add_subdirectory(noise3d)

add_subdirectory(paramgeo2d)
add_subdirectory(paramgeo3d)

add_subdirectory(transform_handler)

add_subdirectory(adf3d)
add_subdirectory(selfcollision_bvh)


add_subdirectory(shader_glsl)
add_subdirectory(offscreen)

add_subdirectory(subdiv)
add_subdirectory(exponential_map)
add_subdirectory(edit_cubes)
add_subdirectory(spherical_harmonics)
add_subdirectory(isosurface_stuffing)
add_subdirectory(mathexpeval)

add_subdirectory(pbd_grid2d)
add_subdirectory(pbd_cloth)

add_subdirectory(cadmsh2d)

add_subdirectory(convhull)
add_subdirectory(dyntet_tetrahedralization)
add_subdirectory(dyntri2d_triangulation)
add_subdirectory(dyntri2d_remesh)
add_subdirectory(dyntri3d_edgecollapse)

add_subdirectory(cad2d)
add_subdirectory(cad3d)

add_subdirectory(fem2d)
add_subdirectory(fem3d)
add_subdirectory(fem25d_eigen)
add_subdirectory(cloth_internal)
add_subdirectory(cloth_selfcollision)

add_subdirectory(dyntri2d_cloth)