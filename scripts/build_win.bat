
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
mkdir ../GTest_Lib
cmake --install . --prefix ../GTest_Lib
cd ../../..

set path_gtest_root=%~dp0..\3rd_party\GTEST_Lib
echo "path_gtest_root: %path_gtest_root%"

cd test_cpp
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 -DUSE_HEADERONLY=ON  -DGTEST_ROOT="%path_gtest_root%"
cmake --build . --config Release
"Release/runUnitTests.exe"
cd ../../

cd test_cpp
mkdir buildVS64Static
cd buildVS64Static
cmake .. -A x64 -DUSE_HEADERONLY=OFF -DGTEST_ROOT="%path_gtest_root%"
cd ../../

: ################################
: build glfw

git submodule update --init 3rd_party/glfw
cd 3rd_party/glfw
cmake . -A x64
cmake --build . --config Release
mkdir ../GLFW_Lib
cmake --install . --prefix ../GLFW_Lib
cd ../..

: ##############################
: glfw_oldgl

set path_glfw_root=%~dp0..\3rd_party\GLFW_Lib
echo "path_glfw_root: %path_glfw_root%"

cd examples_oldgl_glfw
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 -DUSE_HEADERONLY=ON -DGLFW_ROOT=%path_glfw_root%
cmake --build . --config Release
cd ../../

cd examples_oldgl_glfw
mkdir buildVS64Static
cd buildVS64Static
cmake .. -A x64 -DUSE_HEADERONLY=OFF -DGLFW_ROOT=%path_glfw_root%
cmake --build . --config Release
cd ../../

: ##############################
: glfw_newgl

set path_glfw_root=%~dp0..\3rd_party\GLFW_Lib
echo "path_glfw_root: %path_glfw_root%"

cd examples_newgl_glfw
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 -DUSE_HEADERONLY=ON -DGLFW_ROOT=%path_glfw_root%
cmake --build . --config Release
cd ../../

cd examples_newgl_glfw
mkdir buildVS64Static
cd buildVS64Static
cmake .. -A x64 -DUSE_HEADERONLY=OFF -DGLFW_ROOT=%path_glfw_root%
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

set path_glfw_root=%~dp0..\3rd_party\GLFW_Lib
set path_zlib_lib=%~dp0..\3rd_party\zlib\buildMake\Release\zlib.lib
set path_zlib_inc=%~dp0..\3rd_party\zlib
echo "path_glfw_root: %path_glfw_root%"
cd examples_oldgl_glfw_cnpy
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 ^
  -DUSE_HEADERONLY=ON ^
  -DZLIB_LIBRARY="%path_zlib_lib%" ^
  -DZLIB_INCLUDE_DIR="%path_zlib_inc%" ^
  -DGLFW_ROOT=%path_glfw_root%
cmake --build . --config Release
cd ../../

: ##############################
: build OpenEXR

set path_zlib_library=%~dp0..\3rd_party\zlib\buildMake\Release\zlib.lib
set path_zlib_root=%~dp0..\3rd_party\zlib
set path3rdparty=%~dp0..\3rd_party
cd 3rd_party
curl -L https://github.com/AcademySoftwareFoundation/openexr/archive/v2.5.2.zip -o openexr.zip
7z x openexr.zip -y
cd openexr-2.5.2
cmake . -A x64 ^ 
  -DZLIB_ROOT=%path_zlib_root% ^
  -DZLIB_LIBRARY=%path_zlib_library% ^
  -DPYILMBASE_ENABLE=OFF ^ 
  -DBUILD_TESTING=OFF
cmake --build . --config Release
mkdir %path3rdparty%\OpenEXR_Lib
cmake --install . --prefix %path3rdparty%\OpenEXR_Lib
cd ../../

: ##############################
: build Alembic

set path3rdparty=%~dp0..\3rd_party

git submodule update --init 3rd_party/alembic

cd 3rd_party/alembic
git checkout master
git pull origin master
cmake -A x64 . ^
  -DUSE_TESTS=OFF ^
  -DALEMBIC_SHARED_LIBS=ON ^
  -DILMBASE_ROOT=%path3rdparty%\OpenEXR_Lib
cmake --build . --config Release
mkdir %path3rdparty%\Alembic_Lib
cmake --install . --prefix %path3rdparty%\Alembic_Lib
cd ../../

: ###############################
: build examples alembic

cd examples_alembic
mkdir buildVS64Static
cd buildVS64Static
cmake -A x64 .. -DILMBASE_ROOT=%path3rdparty%\OpenEXR_Lib
cmake --build .
cd ../..