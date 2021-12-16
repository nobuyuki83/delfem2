cmake_minimum_required(VERSION 3.12)

####################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
ENDIF()

####################################

project(02_rig)

add_definitions(-DPATH_INPUT_DIR="${PROJECT_SOURCE_DIR}/../../test_inputs")

set(CMAKE_PREFIX_PATH ../../3rd_party/Imathlib)
FIND_PACKAGE(Imath)

set(CMAKE_PREFIX_PATH ../../3rd_party/alembiclib)
FIND_PACKAGE(Alembic)

# delfem2
set(DELFEM2_INCLUDE_DIR "../../include")

####################################

add_executable(${PROJECT_NAME}
	main.cpp
)

target_include_directories(${PROJECT_NAME} 
PUBLIC 
    Alembic::Alembic
    ${DELFEM2_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} 
PUBLIC
    Alembic::Alembic
)


