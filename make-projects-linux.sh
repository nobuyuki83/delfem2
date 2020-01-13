#################################
# download & build submodules

git submodule update --init --recursive

######################
## test

cd 3rd_party/googletest
cmake .
make
cd ../..

cd test_cpp
mkdir buildMake
cd buildMake
cmake ..
make
./runUnitTests
cd ../../

#######################

cd 3rd_party/glfw
cmake .
make
cd ../..

cd examples_glfw
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../

cd examples_glfw_oldgl
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../

#######################
# glut

cd examples_glut
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../

#######################
# cuda

cd examples_cuda
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../

#######################


#virtualenv --python=python3.7 myenv

#PATH_PYTHON="myenv/bin/"
#PATH_PYTHON=$( cd ${PATH_PYTHON}; pwd )"/python3"
PATH_PYTHON=$(which python3)
echo ${PATH_PYTHON}

cd src_pybind/core
mkdir buildMake
cd buildMake
cmake -DPYTHON_EXECUTABLE:PATH=${PATH_PYTHON}  ..
cd ../../../
cd src_pybind/gl

mkdir buildMake
cd buildMake
cmake -DPYTHON_EXECUTABLE:PATH=${PATH_PYTHON}  ..
cd ../../../

pip3 uninstall PyDelFEM2 -y
pip3 uninstall PyDelFEM2 -y
pip3 install -e .
python3 setup.py test

#python3 setup.py install
#python3 setup.py test
#python3 setup.py sdist bdist_wheel
#twine upload dist/*

########################
## build dll for CSharp

#cd src_dll
#mkdir build
#cd build
#cmake ..
#make

########################
