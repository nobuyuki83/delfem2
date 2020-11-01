cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" ON)

IF(NOT USE_HEADERONLY)
	message("$$$$$$$$$$$$$$$$$$$$$")
	message("generate project as STATIC library: USE_HEADERONLY=${USE_HEADERONLY}")
	project(examples_smpl_static)
    add_subdirectory(dfm2_static_glfwold)
ELSE()
	message("$$$$$$$$$$$$$$$$$$$$$")
	message("generate project as HEADER-ONLY library:  USE_HEADERONLY=${USE_HEADERONLY}")
	project(examples_smpl_hdronly)
ENDIF()

add_subdirectory(00_RigSmplArticulated)
add_subdirectory(01_RigSmplTransfer)
add_subdirectory(02_RigSmplIk)
add_subdirectory(03_RigSmplIkArap)
add_subdirectory(04_RigSmplIkImage)
add_subdirectory(05_RigSmplUi)

add_subdirectory(10_RigSmplCloth)
add_subdirectory(11_RigSmplClothPose)
add_subdirectory(12_RigSmplIkImageCloth)


