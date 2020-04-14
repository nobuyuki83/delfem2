cmake_minimum_required(VERSION 2.8)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

######################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

##########################################

project(610_DefArapEdge)

set(DELFEM2_INCLUDE_DIR "../../include")
set(DELFEM2_INC         "../../include/delfem2")
set(3RD_PARTY           "../../3rd_party")
add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

set(GLFW_DIR "${3RD_PARTY}/glfw")
get_filename_component(GLFW_DIR ${GLFW_DIR} ABSOLUTE)

find_package(OpenGL REQUIRED)

#######################################

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLFW_DIR}/include
  ${DELFEM2_INCLUDE_DIR})
  
IF(NOT USE_HEADERONLY)
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  set(DFM2_LIB
      ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/libdfm2_static_glfwold.a )
ENDIF()

add_executable(${PROJECT_NAME}
	main.cpp
)

###############################################

IF(APPLE)
	find_library( COCOA_LIBRARY Cocoa )
	find_library( IOKIT_LIBRARY IOKit )
	find_library( COREVIDEO_LIBRARY CoreVideo )
	target_link_libraries(${PROJECT_NAME}
      ${DFM2_LIB}
			${OPENGL_LIBRARY}
			${GLFW_DIR}/src/libglfw3.a
			${COCOA_LIBRARY}
			${IOKIT_LIBRARY}
			${COREVIDEO_LIBRARY}
			)
ELSEIF(UNIX)
	target_link_libraries(${PROJECT_NAME}
      ${DFM2_LIB}
			${OPENGL_LIBRARY}
			${GLFW_DIR}/src/libglfw3.a
			-lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
			)
ELSEIF(MSVC)
	if(CMAKE_SIZEOF_VOID_P EQUAL 8)
		target_link_libraries(${PROJECT_NAME}
				${GLFW_DIR}/src/Release/glfw3.lib
				${OPENGL_LIBRARY})
	else()
		target_link_libraries(${PROJECT_NAME}
				${GLFW_DIR}/src/Release/glfw3.lib
				${OPENGL_LIBRARY})
	endif()
ENDIF()

