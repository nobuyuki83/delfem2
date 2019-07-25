
git submodule update --init --recursive
cd tests/googletest
cmake .
make
cd ../..

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

cd examples_cpp
mkdir buildXcode
cd buildXcode
cmake -G Xcode ..
cd ../../

cd tests
mkdir buildXcode
cd buildXcode
cmake -G Xcode ..
cd ../../


#python3 setup.py install
#python3 setup.py test
#python3 setup.py sdist bdist_wheel
#twine upload dist/*
