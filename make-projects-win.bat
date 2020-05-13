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

cd examples_glfwold
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake -A x64 -DUSE_HEADERONLY=ON ..
cmake --build . --config Release
cd ../../

cd examples_glfwold
mkdir buildVS64Static
cd buildVS64Static
cmake -A x64 -DUSE_HEADERONLY=OFF ..
cmake --build . --config Release
cd ../../

: ##############################
: glfw_newgl

cd examples_glfwnew
mkdir buildVS64Hdronly
cd buildVS64Hdronly
cmake -A x64 -DUSE_HEADERONLY=ON ..
cmake --build . --config Release
cd ../../

cd examples_glfwnew
mkdir buildVS64Static
cd buildVS64Static
cmake -A x64 -DUSE_HEADERONLY=OFF ..
cmake --build . --config Release
cd ../../

: ###############################
: pybind

cd src_pybind/core
mkdir buildVS64
cd buildVS64
cmake -A x64 -DUSE_HEADERONLY=OFF ..
cmake --build . --config Release
cd ../../../

cd src_pybind/gl
mkdir buildVS64
cd buildVS64
cmake -A x64 -DUSE_HEADERONLY=OFF ..
cmake --build . --config Release
cd ../../../

pip uninstall PyDelFEM2 -y
pip uninstall PyDelFEM2 -y
pip3 install -e .
python setup.py test
