cmake_minimum_required(VERSION 2.8)

project(examples_glut)

# opengl
add_subdirectory(00_viewer3d)
add_subdirectory(01_shader)
add_subdirectory(02_noise2d)
add_subdirectory(03_noise3d)

add_subdirectory(offscreen)
add_subdirectory(glbuffer)

# geometry
add_subdirectory(10_convhull)

add_subdirectory(subdiv)
add_subdirectory(exponentialmap)
add_subdirectory(edit_cubes)
add_subdirectory(spherical_harmonics)
add_subdirectory(meshslice)

add_subdirectory(paramgeo2d)
add_subdirectory(paramgeo3d)
add_subdirectory(transform_handler)

add_subdirectory(selfcollision_bvh)


# meshgen
add_subdirectory(adf3d)
add_subdirectory(iss)

add_subdirectory(dyntet_tetrahedralization)
add_subdirectory(dyntri2d_triangulation)
add_subdirectory(dyntri2d_remesh)
add_subdirectory(dyntri3d_edgecollapse)

# cad
add_subdirectory(cad2d)
add_subdirectory(cad3d)

# fem
add_subdirectory(rigidbody3d_rotation)

add_subdirectory(fem2d)
add_subdirectory(fem2d_helmholtz)
add_subdirectory(fem3d)
add_subdirectory(fem25d_eigen)
add_subdirectory(fem25d_stiffwarp)
add_subdirectory(femcloth_internal)
add_subdirectory(femcloth_selfcollision)
add_subdirectory(femcloth_dtri)
add_subdirectory(fempb_mitc3)
add_subdirectory(fempb_mitc3eigen)

# io
add_subdirectory(bvh_read)
add_subdirectory(gltf_read)

# pbd
add_subdirectory(pbd_grid2d)
add_subdirectory(pbdcloth)
add_subdirectory(pbdcloth_cad)


add_subdirectory(cadmsh2d)
add_subdirectory(pbdcloth_cad_gltf)
