#################################
# download & build submodules

git submodule update --init --recursive

################################
# build examples_glfw

cd 3rd_party/glfw
cmake .
make
cd ../..

cd examples_glfw
mkdir buildXcode
cd buildXcode
cmake -G Xcode ..
cmake --build .
cd ../../

cd examples_glfw
mkdir buildEm
cd buildEm
cmake -DEMSCRIPTEN ..
make
cd ../../


################################
# build examples_glut

cd examples_glfw_oldgl
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../

cd examples_glfw_oldgl
mkdir buildXcode
cd buildXcode
cmake -G Xcode ..
cmake --build .
cd ../../

################################
# build examples_glut

cd examples_glut
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../

cd examples_glut
mkdir buildXcode
cd buildXcode
cmake -G Xcode ..
cmake --build .
cd ../../


##################################
## build CSharp

#cd src_dll 
#mkdir build 
#cd build
#cmake ..
#make
#cd ../../

#cd examples_cs
#mcs helloworld.cs -define:__APPLE__ -out:helloworld.exe
#mono helloworld.exe
#./helloworld.exe
#cd ../

################################
# build python

#virtualenv --python=python3.7 myenv
#PATH_PYTHON="myenv/bin/"
#PATH_PYTHON=$( cd ${PATH_PYTHON}; pwd )"/python3"
PATH_PYTHON=$(which python)
echo ${PATH_PYTHON}

cd src_pybind/core
mkdir buildXcode 
cd buildXcode
cmake -G Xcode -DPYTHON_EXECUTABLE:PATH=${PATH_PYTHON}  ..
cd ../../../

cd src_pybind/gl
mkdir buildXcode 
cd buildXcode
cmake -G Xcode -DPYTHON_EXECUTABLE:PATH=${PATH_PYTHON}  ..
cd ../../../

pip3 uninstall PyDelFEM2 -y
pip3 uninstall PyDelFEM2 -y
pip3 install -e .
#python3 setup.py install
python3 setup.py test
#python3 setup.py sdist bdist_wheel
#twine upload dist/*


################################
# test cpp
# (this takes time so put it in the end)

cd "3rd_party/googletest"
cmake .
make
cd ../../

cd test_cpp
mkdir buildXcode
cd buildXcode
cmake -G Xcode ..
cd ../../

cd test_cpp
mkdir buildMake
cd buildMake
cmake ..
make
./runUnitTests
cd ../../

#################################
# torch extension

pip3 uninstall torch_delfem2 -y
pip3 install torch
cd src_pybind/torch
python3 setup.py develop
cd ../../




