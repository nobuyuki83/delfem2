cmake_minimum_required(VERSION 3.12)

option(USE_STATIC_LIB "compile delfem2 as a static library?" OFF)

####################################

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")
ENDIF ()

####################################

project(41_fem2d_poisson)

# glad
get_filename_component(GLAD_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/../../3rd_party/glad4/include ABSOLUTE)

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
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwnew/dfm2_static_glfwnew.lib)
  else ()
    set(DFM2_LIBRARY
        ${PROJECT_SOURCE_DIR}/../dfm2_static_glfwnew/libdfm2_static_glfwnew.a)
  endif ()
endif ()
set(DELFEM2_INCLUDE_DIR "../../include")

#############################################

include_directories(
    ${GLAD_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
    ../../3rd_party/glad4/include/glad/glad.h
    ../../3rd_party/glad4/src/glad.c
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    OpenGL::GL
    glfw
    )