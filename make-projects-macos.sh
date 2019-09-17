
#################################
# download & build submodules
git submodule update --init --recursive


################################
# test cpp

cd test_cpp/googletest
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

################################
# build examples_glfw

cd examples_glfw/glfw
cmake .
make
cd ../..

cd examples_glfw
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

cd src_pybind/eigen
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
