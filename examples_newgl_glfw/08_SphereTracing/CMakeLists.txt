cmake_minimum_required(VERSION 3.12)

#######################

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
endif ()

############################################

project(08_SphereTracing)

add_definitions(-DPATH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# glad
get_filename_component(GLAD_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/include ABSOLUTE)
SET(GLAD_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/include/glad/glad.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/glad4/src/glad.c)

# glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

# opengl
find_package(OpenGL REQUIRED)

# dfm2
set(DELFEM2_INCLUDE_DIR "../../include")

############################################

include_directories(
    ${GLAD_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
    ${GLAD_SOURCES}
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    OpenGL::GL
    glfw
    )
