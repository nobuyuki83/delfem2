# DelFEM2 C++ Examples using std::thread & Legacy OpenGL

These demos use OpenGL version 2.1 and GLSL shaer version 1.2 which are depricated in many environment. But still it is convenient to use legacy functions such as glBegin(), glEnd(). We will eventually consider porting these demo into newer OpenGL >= 3.3 in the [examples_glfwnew folder](../examples_glfwnew).



## Download dependencies

```bash
# move to the top directory
cd delfem2

# download glfw and compile it
git submodule update --init -- 3rd_party/glfw
cd 3rd_party/glfw
cmake .
make
cd ../..
```





### [10_DefArap](10_DefArap)
<img src="10_DefArap/thumbnail.png" width=200>

