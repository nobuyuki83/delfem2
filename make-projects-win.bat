

cd src_pybind/core
mkdir buildVS
cd buildVS
cmake -A x64 ..
cmake --build . --config Release
cd ../../../


cd src_pybind/gl
mkdir buildVS
cd buildVS
cmake -A x64 -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/glew" ..
cmake --build . --config Release
cd ../../../


cd src_pybind/eigen
mkdir buildVS
cd buildVS
cmake -A x64 -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/glew" ..
cmake --build . --config Release
cd ../../../


cd examples_cpp
mkdir buildVS
cd buildVS
cmake -A Win32 -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/glew;C:/Program Files (x86)/freeglut" ..
cmake --build .
cd ../../

