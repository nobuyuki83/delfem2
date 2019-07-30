

git submodule update --init --recursive
cd src_cpp/external/glfw
mkdir buildVS64
cd buildVS64
cmake -A x64 ..
cmake --build . --config Release
:: cmake --build . --config Debug
cd ..
mkdir buildVS32
cd buildVS32
cmake -A Win32 ..
cmake --build . --config Release
:: cmake --build . --config Debug
cd ../
cd ../../../





cd src_pybind/core
mkdir buildVS64
cd buildVS64
cmake -A x64 ..
cmake --build . --config Release
cd ../../../

cd src_pybind/gl
mkdir buildVS64
cd buildVS64
cmake -A x64 -DCMAKE_PREFIX_PATH="C:/Program Files/glew" ..
cmake --build . --config Release
cd ../../../


cd src_pybind/eigen
mkdir buildVS64
cd buildVS64
cmake -A x64 -DCMAKE_PREFIX_PATH="C:/Program Files/glew" ..
cmake --build . --config Release
cd ../../../


cd examples_cpp
mkdir buildVS32
cd buildVS32
cmake -A Win32 -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/glew;C:/Program Files (x86)/freeglut" ..
cmake --build .
cd ../../


cd examples_cpp
mkdir buildVS64
cd buildVS64
cmake -A x64 -DCMAKE_PREFIX_PATH="C:/Program Files/glew;C:/Program Files/freeglut" ..
cmake --build .
cd ../../



pip uninstall PyDelFEM2 -y
pip uninstall PyDelFEM2 -y
python setup.py install
python setup.py test



cd tests\googletest
mkdir buildVS64
cd buildVS64
cmake -A x64 .. 
cmake --build . --config Release
cd ..\..\..
cd tests
mkdir buildVS64
cd buildVS64
cmake -A x64 ..
cmake --build . --config Release
Release\runUnitTests.exe
cd ..\..\