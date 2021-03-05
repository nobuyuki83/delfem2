
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

git submodule update --init 3rd_party/glfw
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

git submodule update --init 3rd_party/zlib
cd 3rd_party/zlib
mkdir buildMake
cd buildMake
cmake -A x64 ..
cmake --build . --config Release
copy zconf.h ..
cd ../../../

: ##############################
: glfw_smpl

cd examples_oldgl_glfw_cnpy
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake .. -A x64 ^
  -DUSE_HEADERONLY=ON ^
  -DZLIB_LIBRARY=%path3rdparty%\zlib\buildMake\Release\zlib.lib ^
  -DZLIB_INCLUDE_DIR="%path3rdparty%\zlib"
cmake --build . --config Release
cd ../../

: ##############################
: build OpenEXR

cd 3rd_party
curl -L https://github.com/AcademySoftwareFoundation/openexr/archive/v2.5.2.zip -o openexr.zip
7z x openexr.zip -y
cd openexr-2.5.2
cmake -A x64 . -DZLIB_ROOT=%zlib_root% -DZLIB_LIBRARY=%zlib_library% -DPYILMBASE_ENABLE=OFF -DBUILD_TESTING=OFF
cmake --build . --config Release
mkdir %path3rdparty%\OpenEXR
cmake --install . --prefix %path3rdparty%\OpenEXR
cd ../../

: ##############################
: build Alembic

git submodule update --init 3rd_party/alembic
cd 3rd_party/alembic
git checkout master
git pull origin master
: cmake -A x64 . -DUSE_TESTS=OFF -DALEMBIC_SHARED_LIBS=ON -DALEMBIC_ILMBASE_ROOT=%pathopenexr%\bin\Release -DILMBASE_INCLUDE_DIR=%pathopenexr%\bin\Release -DILMBASE_ROOT=%pathopenexr%\bin\Release
: cmake -A x64 . -DUSE_TESTS=OFF -DALEMBIC_SHARED_LIBS=ON -DALEMBIC_ILMBASE_ROOT=%pathopenexr%\bin\Release -DILMBASE_ROOT=%pathopenexr%\bin\Release
cmake -A x64 . ^
  -DUSE_TESTS=OFF ^
  -DALEMBIC_SHARED_LIBS=ON ^
  -DILMBASE_ROOT=%path3rdparty%\OpenEXR 
cmake --build . --config Release
cd ../../