[![wercker status](https://app.wercker.com/status/03b6d924ec82270e22a04c3584fbf4de/s/master "wercker status")](https://app.wercker.com/project/byKey/03b6d924ec82270e22a04c3584fbf4de)  [![travis_status](https://travis-ci.org/nobuyuki83/delfem2.svg?branch=master)](https://travis-ci.org/nobuyuki83/delfem2)


# DelFEM2

A handy toolset for coding geometry processing and fem simulation

The implementation is based on the [DelFEM](https://github.com/nobuyuki83/DelFEM) library

Please find out more detail in this [project document](https://nobuyuki83.github.io/delfem2/)


# Dependency

PyDelFEM is run on Python3. Python2 is not supported.

PyDelFEM2 depends on following python packages:
- Numpy
- glfw
- PyOpenGL  
- PySide2
These packages are automatically installed when you installed the PyDelFEM2 using the ```setup.py``` or ```pip3```.

# Install

## from PyPl

PyDelFEM2 can be installed simply with 

```
pip3 install PyDelFEM2
```

If you don't have pip install it with

```
sudo apt-get install python3-pip
```


Installation may fails because the necessary OpenGL packages are missing. Run following commanad for Ubuntu.

```
sudo apt-get install freeglut3-dev libglfw3-dev libglew-dev
```



## from GitHub

installation from github can be done with:
```
pip3 install git+https://github.com/nobuyuki83/delfem2
```


## building from the source code

```
git clone https://github.com/nobuyuki83/delfem2.git
git submodle update --init --recursive
python3 setup.py install
```

if you don't have git, install it with 
```
sudo apt-get install git
```



