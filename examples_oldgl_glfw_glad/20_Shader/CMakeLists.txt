cmake_minimum_required(VERSION 3.12)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
endif ()

#######################

project(20_Shader)

add_definitions(-DPATH_ASSET_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")
add_definitions(-DPATH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# glad
get_filename_component(GLAD_INCLUDE_DIR ../../3rd_party/glad3/include ABSOLUTE)

# glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

# opengl
find_package(OpenGL REQUIRED)

# delfem2
set(DELFEM2_INCLUDE_DIR "../../include")

#####################################

include_directories(
    ${GLAD_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party
)

set(3RD_PARTY "../../3rd_party")
add_executable(${PROJECT_NAME}
    ${3RD_PARTY}/glad3/include/glad/glad.h
    ${3RD_PARTY}/glad3/src/glad.c
    main.cpp
    )

#####################################

target_link_libraries(${PROJECT_NAME}
    OpenGL::GL
    glfw
    )