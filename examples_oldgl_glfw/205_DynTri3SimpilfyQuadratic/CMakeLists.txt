cmake_minimum_required(VERSION 3.12)

option(USE_STATIC_LIB "compile delfem2 with static library?" OFF)

######################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "/W4 /EHsc")
ELSE ()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
ENDIF ()

####################################

project(205_DynTri3SimplifyQuadratic)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
IF (USE_STATIC_LIB)
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  if (MSVC)
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/dfm2_static_glfwold.lib)
  else ()
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/libdfm2_static_glfwold.a)
  endif ()
ENDIF ()
set(DELFEM2_INCLUDE_DIR "../../include")

#########################################

include_directories(        
    ${DELFEM2_INCLUDE_DIR}
    )

add_executable(${PROJECT_NAME}
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    OpenGL::GL
    glfw
    )
