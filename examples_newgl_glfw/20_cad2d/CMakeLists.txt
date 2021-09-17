cmake_minimum_required(VERSION 3.12)

option(USE_STATIC_LIB "compile delfem2 as a static library?" OFF)

###########################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
ENDIF ()


###########################

project(20_cad2d)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

###########################

# glfw
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}../../3rd_party/libglfw)
find_package(glfw3 REQUIRED)

# opengl
find_package(OpenGL REQUIRED)

# glad
set(GLAD_INCLUDE_DIR "../../3rd_party/glad4/include")
get_filename_component(GLAD_INCLUDE_DIR ${GLAD_INCLUDE_DIR} ABSOLUTE)

# dfm2
IF (USE_STATIC_LIB)
  add_definitions(-DDFM2_STATIC_LIBRARY=ON)
  if (MSVC)
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwnew/dfm2_static_glfwnew.lib)
  else ()
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwnew/libdfm2_static_glfwnew.a)
  endif ()
endif ()
set(DELFEM2_INCLUDE_DIR "../../include")

########################

include_directories(
    ${GLAD_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

set(3RD_PARTY "../../3rd_party")
add_executable(${PROJECT_NAME}
    ${3RD_PARTY}/glad4/include/glad/glad.h
    ${3RD_PARTY}/glad4/src/glad.c
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    OpenGL::GL
    glfw
    )
