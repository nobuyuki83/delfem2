cmake_minimum_required(VERSION 3.12)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

IF(NOT USE_HEADERONLY)
	message("$$$$$$$$$$$$$$$$$$$$$")
	message("generate project as STATIC library: USE_HEADERONLY=${USE_HEADERONLY}")
	project(examples_oldgl_glfw_cnpy_static)
    add_subdirectory(dfm2_static_glfwold)
ELSE()
	message("$$$$$$$$$$$$$$$$$$$$$")
	message("generate project as HEADER-ONLY library:  USE_HEADERONLY=${USE_HEADERONLY}")
	project(examples_oldgl_glfw_cnpy_hdronly)
ENDIF()

add_subdirectory(00_PoseBone)
add_subdirectory(01_RigTransfer)
add_subdirectory(02_Ik)
add_subdirectory(03_IkArap)
add_subdirectory(04_IkImage)
add_subdirectory(05_Ui)

add_subdirectory(10_Cloth)
add_subdirectory(11_ClothPose)
add_subdirectory(12_IkImageCloth)
add_subdirectory(13_ClothPoseRig)
add_subdirectory(14_ClothPoseTexture)
add_subdirectory(15_ClothInitialConfig)

add_subdirectory(20_PoseBoneBlendshape)

