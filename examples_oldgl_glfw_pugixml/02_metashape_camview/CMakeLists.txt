cmake_minimum_required(VERSION 3.12)

######################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
ENDIF ()

######################################################

project(02_metashape_camview)

add_definitions(-DPATH_SOURCE_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

# glfw
set(CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

# opengl
find_package(OpenGL REQUIRED)

set(DELFEM2_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../include)

######################################################

include_directories(
    ${DELFEM2_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/../../3rd_party/pugixml/src
    ${PROJECT_SOURCE_DIR}/../../3rd_party
)

add_executable(${PROJECT_NAME}
    main.cpp
    ${PROJECT_SOURCE_DIR}/../../3rd_party/pugixml/src/pugixml.cpp
)

target_link_libraries(${PROJECT_NAME}
    OpenGL::GL
    glfw
    )