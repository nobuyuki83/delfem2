
git submodule update --init --recursive

: ################################
: test cpp

git submodule update --init -- 3rd_party/googletest
cd 3rd_party/googletest
git checkout master
git pull origin master
mkdir buildVS64
cd buildVS64
cmake -A x64 -Dgtest_force_shared_crt=ON ..
cmake --build . --config Release
cmake --install . --prefix ../../libgtest
cd ../../..

set path_gtest_root=%~dp0..\3rd_party\libgtest
echo "path_gtest_root: %path_gtest_root%"

cd test_cpp
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 -DUSE_STATIC_LIB=OFF  -DGTEST_ROOT="%path_gtest_root%"
cmake --build . --config Release
"Release/runUnitTests.exe"
cd ../../

cd test_cpp
mkdir buildVS64Static
cd buildVS64Static
cmake .. -A x64 -DUSE_STATIC_LIB=ON -DGTEST_ROOT="%path_gtest_root%"
cmake --build . --config Release
cd ../../

: ################################
: build glfw

git submodule update --init 3rd_party/glfw
cd 3rd_party/glfw
cmake . -A x64
cmake --build . --config Release
cmake --install . --prefix ../libglfw
cd ../..

: ##############################
: glfw_oldgl

set path_glfw_root=%~dp0..\3rd_party\libglfw
echo "path_glfw_root: %path_glfw_root%"

cd examples_oldgl_glfw
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 -DUSE_STATIC_LIB=OFF
cmake --build . --config Release
cd ../../

cd examples_oldgl_glfw
mkdir buildVS64Static
cd buildVS64Static
cmake .. -A x64 -DUSE_STATIC_LIB=ON
cmake --build . --config Release
cd ../../

: ##############################
: glfw_newgl

set path_glfw_root=%~dp0..\3rd_party\libglfw
echo "path_glfw_root: %path_glfw_root%"

cd examples_newgl_glfw
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 -DUSE_STATIC_LIB=OFF
cmake --build . --config Release
cd ../../

cd examples_newgl_glfw
mkdir buildVS64Static
cd buildVS64Static
cmake .. -A x64 -DUSE_STATIC_LIB=ON
cmake --build . --config Release
cd ../../

: ##############################
: build zlib

git submodule update --init 3rd_party/zlib
cd 3rd_party/zlib
mkdir buildMake
cd buildMake
cmake -A x64 ..
cmake --build . --config Release
copy zconf.h ..
cd ../../../

: ##############################
: glfw_cnpy

set path_glfw_root=%~dp0..\3rd_party\libglfw
set path_zlib_lib=%~dp0..\3rd_party\zlib\buildMake\Release\zlib.lib
set path_zlib_inc=%~dp0..\3rd_party\zlib
echo "path_glfw_root: %path_glfw_root%"
cd examples_oldgl_glfw_cnpy
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 ^
  -DZLIB_LIBRARY="%path_zlib_lib%" ^
  -DZLIB_INCLUDE_DIR="%path_zlib_inc%" ^
  -DGLFW_ROOT=%path_glfw_root%
cmake --build . --config Release
cd ../../

: ##############################
: build Imath

git submodule update --init -- 3rd_party/Imath
cd 3rd_party/Imath 
git checkout main
git pull origin main
mkdir build
cd build 
cmake -A x64 .. ^
  -DBUILD_SHARED_LIBS=OFF ^
  -DBUILD_TESTING=OFF
cmake --build . --config Release
cmake --install . --prefix ../../Imathlib
cd ../../../

: ##############################
: build Alembic

set path3rdparty=%~dp0..\3rd_party

git submodule update --init 3rd_party/alembic

cd 3rd_party/alembic
git checkout master
git pull origin master
mkdir build
cd build
cmake -A x64 .. ^
  -DUSE_TESTS=OFF ^
  -DALEMBIC_SHARED_LIBS=OFF ^
  -DCMAKE_PREFIX_PATH="%path3rdparty%\Imathlib"
cmake --build . --config Release
cmake --install . --prefix ../../alembiclib
cd ../../../

: ###############################
: build examples alembic

cd examples_alembic
mkdir buildVS64
cd buildVS64
cmake -A x64 ..
cmake --build . --config Release
cd ../..