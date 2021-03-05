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
copy zconf.h ..
cd ../../../
set zlib_root=%~dp0..\3rd_party\zlib\
set zlib_library=%zlib_root%\buildMake\Release\zlib.lib
echo "zlib_root: %zlib_root%"
echo "zlib_library: %zlib_library%"

: ##############################
: glfw_smpl

cd examples_oldgl_glfw_cnpy
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake -A x64 -DUSE_HEADERONLY=ON -DZLIB_LIBRARY=%zlib_library% -DZLIB_INCLUDE_DIR="..\..\3rd_party\zlib;..\..\3rd_party\zlib\buildMake" ..
cmake --build . --config Release
cd ../../



: ##############################

goto :END
: ##############################
: build OpenEXR

cd 3rd_party
curl -L https://github.com/AcademySoftwareFoundation/openexr/archive/v2.5.2.zip -o openexr.zip
7z x openexr.zip -y
cd openexr-2.5.2
cmake -A x64 . -DZLIB_ROOT=%zlib_root% -DZLIB_LIBRARY=%zlib_library% -DPYILMBASE_ENABLE=OFF -DBUILD_TESTING=OFF
cmake --build . --config Release
cd ../../
set pathopenexr=%~dp0..\3rd_party\openexr-2.5.2
copy %pathopenexr%\IlmBase\Half\Release\*.lib %pathopenexr%\bin\Release
copy %pathopenexr%\IlmBase\Iex\Release\*.lib %pathopenexr%\bin\Release
copy %pathopenexr%\IlmBase\IexMath\Release\*.lib %pathopenexr%\bin\Release
copy %pathopenexr%\IlmBase\IlmThread\Release\*.lib %pathopenexr%\bin\Release
copy %pathopenexr%\IlmBase\IMath\Release\*.lib %pathopenexr%\bin\Release
copy %pathopenexr%\IlmBase\config\IlmBaseConfig.h %pathopenexr%\bin\Release

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
  -DILMBASE_ROOT=%pathopenexr%\bin\Release ^
  -DILMBASE_INCLUDE_DIR=%pathopenexr%\IlmBase\Half
:  -DALEMBIC_ILMBASE_ROOT=%pathopenexr%\bin\Release
cmake --build . --config Release
cd ../../

exit /b
:END