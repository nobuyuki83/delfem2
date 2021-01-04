git submodule update --init --recursive

: ################################
: test cpp

cd 3rd_party/googletest
mkdir buildVS64
cd buildVS64
cmake -A x64 -Dgtest_force_shared_crt=ON ..
cmake --build . --config Release
cd ../../..

cd test_cpp
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake -A x64 -DUSE_HEADERONLY=ON ..
cmake --build . --config Release
"Release/runUnitTests.exe"
cd ../../

cd test_cpp
mkdir buildVS64Static
cd buildVS64Static
cmake -A x64 -DUSE_HEADERONLY=OFF ..
cd ../../

: ################################
: build glfw

cd 3rd_party/glfw
cmake . -A x64
cmake --build . --config Release
cd ../..

: ##############################
: glfw_oldgl

cd examples_oldgl_glfw
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake -A x64 -DUSE_HEADERONLY=ON ..
cmake --build . --config Release
cd ../../

cd examples_oldgl_glfw
mkdir buildVS64Static
cd buildVS64Static
cmake -A x64 -DUSE_HEADERONLY=OFF ..
cmake --build . --config Release
cd ../../

: ##############################
: glfw_newgl

cd examples_newgl_glfw
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake -A x64 -DUSE_HEADERONLY=ON ..
cmake --build . --config Release
cd ../../

cd examples_newgl_glfw
mkdir buildVS64Static
cd buildVS64Static
cmake -A x64 -DUSE_HEADERONLY=OFF ..
cmake --build . --config Release
cd ../../

: ##############################
: build zlib

cd 3rd_party
git clone https://github.com/madler/zlib.git
cd zlib
mkdir buildMake
cd buildMake
cmake -A x64 ..
cmake --build . --config Release
cd ../../../
set zlib_library=%~dp0..\3rd_party\zlib\buildMake\Release\zlib.lib
echo %zlib_library%

: ##############################
: glfw_smpl

cd examples_oldgl_glfw_cnpy
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake -A x64 -DUSE_HEADERONLY=ON -DZLIB_LIBRARY=%zlib_library% -DZLIB_INCLUDE_DIR="..\..\3rd_party\zlib;..\..\3rd_party\zlib\buildMake" ..
cmake --build . --config Release
cd ../../
