cmake_minimum_required(VERSION 3.12)

option(DFM2_HEADER_ONLY "compile delfem2 as a header-only library?" ON)

#################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if (MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
else ()
  set(CMAKE_CXX_FLAGS "-Wall  -g")
endif ()

#################################

project(632_PbdSpring2)

# opengl
find_package(OpenGL REQUIRED)

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# dfm2
IF (DFM2_HEADER_ONLY)
  add_definitions(-DDFM2_HEADER_ONLY=ON)
else ()
  if (MSVC)
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/dfm2_static_glfwold.lib)
  else ()
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwold/libdfm2_static_glfwold.a)
  endif ()
ENDIF ()
set(DELFEM2_INCLUDE_DIR "../../include")

#################################

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARIES}
    )
