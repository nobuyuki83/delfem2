cmake_minimum_required(VERSION 3.12)

option(DFM2_HEADER_ONLY "compile delfem2 as a header-only library?" ON)

#######################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if (MSVC)
  set(CMAKE_CXX_FLAGS "-O2 /W4 \
     /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4800 /wd4996 \
     /wd4530 /wd4245 /wd4505 /wd4838 /wd4305 /wd4505 /wd4456")
else ()
  set(CMAKE_CXX_FLAGS "-g \
    -Wall -Wno-deprecated-declarations")
endif ()

#######################

project(21_ShaderTexLaplace)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glad
get_filename_component(GLAD_INCLUDE_DIR ../../3rd_party/glad3/include ABSOLUTE)

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/GLFW_Lib)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
if (DFM2_HEADER_ONLY)
  add_definitions(-DDFM2_HEADER_ONLY=ON)
else ()
  if (MSVC)
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/dfm2_static_glfwold.lib)
  else ()
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/libdfm2_static_glfwold.a)
  endif ()
endif ()
set(DELFEM2_INCLUDE_DIR "../../include")


#####################################

set(3RD_PARTY "../../3rd_party")
include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLAD_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
    ${3RD_PARTY}
)

add_executable(${PROJECT_NAME}
    ${3RD_PARTY}/glad3/include/glad/glad.h
    ${3RD_PARTY}/glad3/src/glad.c
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARIES}
    )