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

add_subdirectory(000_RigSmplArticulated)
add_subdirectory(001_RigSmplTransfer)
add_subdirectory(002_RigSmplIk)
add_subdirectory(003_RigSmplIkArap)
add_subdirectory(004_RigSmplIkImage)
add_subdirectory(005_RigSmplIkImageCloth)
add_subdirectory(006_RigSmplUi)
add_subdirectory(007_RigSmplCloth)
add_subdirectory(008_RigSmplClothPose)


