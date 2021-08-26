cmake_minimum_required(VERSION 3.12)

####################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -g")
ENDIF()

####################################

project(01_mpmsnow)

# dfm2
add_definitions(-DDFM2_HEADER_ONLY=ON)
set(DELFEM2_INCLUDE_DIR "../../include")

#######################################################

include_directories(
    ${DELFEM2_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  main.cpp
)
