cmake_minimum_required(VERSION 3.12)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
ENDIF ()

####################################

project(03_IkArap)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

# opengl
find_package(OpenGL REQUIRED)

# zlib
find_package(ZLIB REQUIRED)

# dfm2
set(DELFEM2_INCLUDE_DIR "../../include")

#############################################

set(3RD_PARTY "../../3rd_party")
include_directories(
    ${DELFEM2_INCLUDE_DIR}
    ${3RD_PARTY}
)

set(DELFEM2_INC "../../include/delfem2")
add_executable(${PROJECT_NAME}
    ${3RD_PARTY}/cnpy/cnpy.h
    ${3RD_PARTY}/cnpy/cnpy.cpp
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    OpenGL::GL
    glfw
    ZLIB::ZLIB
    )