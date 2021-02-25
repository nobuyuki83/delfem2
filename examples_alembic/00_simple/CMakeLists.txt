cmake_minimum_required(VERSION 3.12)

####################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)
IF(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4100 /wd4458 /wd4577 /wd4267 /wd4244 /wd4505 /wd4838 \
    /wd4800 /wd4996 /wd4530 /wd4245 /wd4505 /wd4505 /wd4456 ")
ELSE()
  set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g")
ENDIF()

####################################

project(00_simple)

include(../../3rd_party/alembic/cmake/Modules/FindIlmBase.cmake)
get_filename_component(ALEMBIC_LIBRARY ../../3rd_party/alembic/lib/Alembic/libAlembic.a ABSOLUTE)


####################################

add_executable(${PROJECT_NAME}
	main.cpp
)
target_include_directories(${PROJECT_NAME} 
PUBLIC 
	../../3rd_party/alembic/lib
	${ILMBASE_INCLUDE_DIR}
)
target_link_libraries(${PROJECT_NAME} 
PUBLIC
	${ALEMBIC_LIBRARY}
	${ALEMBIC_ILMBASE_HALF_LIB}
	)


