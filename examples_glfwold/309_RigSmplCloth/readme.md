# 309_PbdClothCadSmpl
![](thumbnail.png)


## How to run the demo:

- download from GitHub and initialize
```
git clone https://github.com/nobuyuki83/delfem2.git
cd delfem2
git submodule update --init --recursive
cd 3rd_party/glfw
cmake .
make
cd ../../
```

- compile the demo
```
mkdir build && cd build
cmake ..
make
```

- get SMPL model
 - get ```pose.txt``` from ```generator_cloth_deformation/img2pose/```
 - put it in the ```generator_cloth_deformation/clothsim/```
 - Then the mesh of clothing will be generated.

