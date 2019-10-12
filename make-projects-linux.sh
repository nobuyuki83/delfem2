
git submodule update --init --recursive

cd test_cpp/googletest
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../..

cd src_cpp/external/glfw
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../../../

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

cd src_pybind/eigen
mkdir buildMake
cd buildMake
cmake -DPYTHON_EXECUTABLE:PATH=${PATH_PYTHON}  ..
cd ../../../

cd tests
mkdir buildMake
cd buildMake
cmake ..
make
./runUnitTests
cd ../../

cd examples_cpp
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../




#python3 setup.py install
#python3 setup.py test
#python3 setup.py sdist bdist_wheel
#twine upload dist/*
