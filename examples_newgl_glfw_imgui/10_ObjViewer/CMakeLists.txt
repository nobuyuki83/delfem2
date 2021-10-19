cmake_minimum_required(VERSION 2.8.12)

#####################################

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(NOT EMSCRIPTEN)
  IF(MSVC)
    set(CMAKE_CXX_FLAGS "/W4")
  ELSE()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
  ENDIF()
ENDIF()

#####################################

project(1_ObjViewer)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

#####################################

# opengl
find_package(OpenGL REQUIRED)

# glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

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
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${GLAD_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party
)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/include/glad/glad.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/src/glad.c
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/ImGuiFileDialog/ImGuiFileDialog.cpp
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    OpenGL::GL
    glfw
    )