cmake_minimum_required(VERSION 2.8)

project(examplesDelFEM2)

add_subdirectory(viewer3d_glut)
add_subdirectory(viewer3d_glfw)

add_subdirectory(noise2d)
add_subdirectory(noise3d)

add_subdirectory(paramgeo2d)
add_subdirectory(paramgeo3d)

add_subdirectory(transform_handler)

add_subdirectory(mathexpeval)

add_subdirectory(bvh_read)
add_subdirectory(gltf_read)

add_subdirectory(shader_glsl)
add_subdirectory(offscreen)
add_subdirectory(glbuffer)

add_subdirectory(subdiv)
add_subdirectory(exponentialmap)
add_subdirectory(edit_cubes)
add_subdirectory(spherical_harmonics)
add_subdirectory(meshslice)
add_subdirectory(convhull)

add_subdirectory(adf3d)
add_subdirectory(iss)

add_subdirectory(dyntet_tetrahedralization)
add_subdirectory(dyntri2d_triangulation)
add_subdirectory(dyntri2d_remesh)
add_subdirectory(dyntri3d_edgecollapse)

add_subdirectory(cad2d)
add_subdirectory(cad3d)

add_subdirectory(fem2d)
add_subdirectory(fem3d)
add_subdirectory(fem25d_eigen)
add_subdirectory(femcloth_internal)
add_subdirectory(femcloth_selfcollision)

add_subdirectory(cadmsh2d)
add_subdirectory(dyntri2d_cloth)

add_subdirectory(selfcollision_bvh)

add_subdirectory(pbd_grid2d)
add_subdirectory(pbdcloth)
add_subdirectory(pbdcloth_cad)
add_subdirectory(pbdcloth_cad_gltf)

add_subdirectory(rigidbody3d_rotation)