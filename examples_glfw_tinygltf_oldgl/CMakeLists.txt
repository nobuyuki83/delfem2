cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

IF(NOT USE_HEADERONLY)
	message("$$$$$$$$$$$$$$$$$$$$$")
	message("generate project as STATIC library: USE_HEADERONLY=${USE_HEADERONLY}")
	project(examples_glfwold_static)
    add_subdirectory(dfm2_static_glfwold)
ELSE()
	message("$$$$$$$$$$$$$$$$$$$$$")
	message("generate project as HEADER-ONLY library:  USE_HEADERONLY=${USE_HEADERONLY}")
	project(examples_glfwold_hdronly)
ENDIF()

# basic funcs
add_subdirectory(00_RigReadFileGLTF)
add_subdirectory(10_PbdClothCadGltf)
