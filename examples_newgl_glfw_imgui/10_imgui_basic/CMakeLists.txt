cmake_minimum_required(VERSION 2.8.12)

#####################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(NOT EMSCRIPTEN)
  IF(MSVC)
    set(CMAKE_CXX_FLAGS "/W4")
  ELSE()
    set(CMAKE_CXX_FLAGS "-Wall -g")
  ENDIF()
ENDIF()

#####################################

project(10_imgui_basic)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

#####################################

# opengl
find_package(OpenGL REQUIRED)

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/FindGLFW.cmake)

# imgui
add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_GLAD)
set(IMGUI_DIR "../../3rd_party/imgui")

# glad
set(GLAD_INCLUDE_DIR "../../3rd_party/glad4/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

# dfm2
add_definitions(-DDFM2_HEADER_ONLY=ON)
set(DELFEM2_INCLUDE_DIR "../../include")

##############################

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${GLAD_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

set(3RD_PARTY           "../../3rd_party")
add_executable(${PROJECT_NAME}
    ${IMGUI_DIR}/imgui.h
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.h
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.h
    ${3RD_PARTY}/glad4/include/glad/glad.h
    ${3RD_PARTY}/glad4/src/glad.c
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
  ${DFM2_LIBRARY}
  ${OPENGL_LIBRARY}
  ${GLFW_LIBRARIES}
)
