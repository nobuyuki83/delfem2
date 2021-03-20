# How to Set Up GLFW Library



Namely, you can set up GLFW in three ways  

- download `glfw` using package manager (for Mac and Ubuntu)
- download pre-build library (Mac and Windows)
- compile the GLFW source code by yourself

Below, we discuss these options in detail 

----



## From Package Manager

for Mac, install `glfw` using package manager `brew` as

```bash
$ brew install glfw
```
For ubuntu, install `glfw` using `apt-get` as

```bash
$ sudo apt-get install -y libx11-dev xorg-dev \
                          libglu1-mesa libglu1-mesa-dev \
                          libgl1-mesa-glx libgl1-mesa-dev
$ sudo apt install -y libglfw3 libglfw3-dev
$ sudo apt install -y libglew-dev
```
Unfortunately, for windows, there is not a easyway to install `glfw` with commands.

---



## Pre-Build Binary Library

You can download the prebuild library from here

 https://www.glfw.org/download.html

Extract the compressed file and rename it as `GLFW_Lib` and put it under the `3rd_party/` folder of the reository. 

Make sure you have a header file `glfw3.h` at

```
ProjectTemplate_CppGlfw/3rd_party/GLFW_Lib/include/GLFW/glfw3.h
```





---


## Build from Source Code

Alternatively, you can build `glfw` from source code and put the library under `3rd_party/GLFW_Lib` with

```bash
$ mkdir 3rd_party/GLFW_Lib 
$ git submodule update --init 3rd_party/glfw
$ cd 3rd_party/glfw
$ cmake .
$ cmake --build . --config Release
$ cmake --install . --prefix ../GLFW_Lib 
```

Make sure you have a header file `glfw3.h` at

```
ProjectTemplate_CppGlfw/3rd_party/GLFW_Lib/include/GLFW/glfw3.h
```







 



