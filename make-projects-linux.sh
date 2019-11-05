#################################
# download & build submodules

git submodule update --init --recursive

######################
## test

cd test_cpp/googletest
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

cd src_cpp/external/glfw
cmake .
make
cd ../../..

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

cd examples_glut
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../


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


#python3 setup.py install
#python3 setup.py test
#python3 setup.py sdist bdist_wheel
#twine upload dist/*
