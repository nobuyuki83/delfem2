
echo "################################"
echo "build examples_glut"
echo "################################"

cd examples_oldgl_glut
mkdir buildMake 
cd buildMake
cmake ..
make
cd ../../

cd examples_oldgl_glut
mkdir buildXcode 
cd buildXcode
cmake -G Xcode ..
cmake --build .
cd ../../


echo "################################"
echo "fetch latest glfw and compile it"
echo "################################"

git submodule update --init -- 3rd_party/glfw
cd 3rd_party/glfw
git checkout master
git pull origin master
cmake .
make
cd ../..


echo "################################"
echo "build examples_glfwold"
echo "################################"

cd examples_oldgl_glfw
mkdir buildXcodeHdronly 
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
# cmake --build . # skip build to save time
cd ../../

cd examples_oldgl_glfw
mkdir buildXcodeStatic 
cd buildXcodeStatic
cmake -G Xcode -DUSE_HEADERONLY=OFF ..
# cmake --build . # skip build to save time
cd ../../

cd examples_oldgl_glfw
mkdir buildMakeHdronly
cd buildMakeHdronly
cmake -DUSE_HEADERONLY=ON ..
make
cd ../../

cd examples_oldgl_glfw
mkdir buildMakeStatic 
cd buildMakeStatic
cmake -DUSE_HEADERONLY=OFF ..
make
cd ../../

echo "################################"
echo "build examples_glfwnew"
echo "################################"

cd examples_newgl_glfw
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
cmake --build .
cd ../../

cd examples_newgl_glfw
mkdir buildXcodeStatic
cd buildXcodeStatic
cmake -G Xcode -DUSE_HEADERONLY=OFF ..
cmake --build .
cd ../../

echo "################################"
echo "build examples_glfw_thread_oldgl"
echo "################################"

cd examples_oldgl_glfw_thread
mkdir buildXcodeHdronly 
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
# cmake --build . # skip build to save time
cd ../../


echo "################################"
echo "fetch latest imgui"
echo "################################"

git submodule update --init -- 3rd_party/imgui
cd 3rd_party/imgui
git checkout master
git pull origin master
cd ../../

echo "################################"
echo "build examples_newgl_glfw_imgui"
echo "################################"

cd examples_newgl_glfw_imgui
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
cmake --build .
cd ../../


echo "################################"
echo "fetch latest tinygltf"
echo "################################"

git submodule update --init -- 3rd_party/tinygltf
cd 3rd_party/tinygltf
git checkout master
git pull origin master
cd ../../

echo "################################"
echo "compile demos using tinygltf"
echo "################################"

cd examples_oldgl_glfw_tinygltf
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
cmake --build .
cd ../../

cd examples_oldgl_glfw_tinygltf
mkdir buildXcodeStatic
cd buildXcodeStatic
cmake -G Xcode -DUSE_HEADERONLY=OFF ..
cmake --build .
cd ../../


echo "################################"
echo "fetch latest cnpy"
echo "################################"

git submodule update --init -- 3rd_party/cnpy
cd 3rd_party/cnpy
git checkout master
git pull origin master
cd ../../

echo "################################"
echo "build examples_cnpy"
echo "################################"

cd examples_oldgl_glfw_cnpy
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
# cmake --build . # skip build to save time
cd ../../

cd examples_oldgl_glfw_cnpy
mkdir buildMakeHdronly
cd buildMakeHdronly
cmake -DUSE_HEADERONLY=ON ..
make 
cd ../../

echo "################################"
echo "build examples_eigen"
echo "################################"

cd examples_oldgl_glfw_eigen
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
cd ../../

cd examples_oldgl_glfw_eigen
mkdir buildMakeHdronly
cd buildMakeHdronly
cmake -DUSE_HEADERONLY=ON ..
make
cd ../../


echo "################################"
echo "build examples_thread"
echo "################################"

cd examples_oldgl_glfw_thread
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
cmake --build .
cd ../../

cd examples_oldgl_glfw_thread
mkdir buildMakeHdronly
cd buildMakeHdronly
cmake -DUSE_HEADERONLY=ON ..
make
cd ../../


echo "######################################"
echo "fetch latest googletest and compile it"
echo "######################################"

git submodule update --init -- 3rd_party/googletest
cd 3rd_party/googletest
git checkout master
git pull origin master
cmake .
make
cd ../../

echo "###############################"
echo "test cpp"
echo "###############################"

cd test_cpp
mkdir buildXcodeStatic
cd buildXcodeStatic
cmake -G Xcode -DUSE_HEADERONLY=OFF ..
cd ../../

cd test_cpp
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
cd ../../

cd test_cpp
mkdir buildMakeHdronly
cd buildMakeHdronly
cmake -DUSE_HEADERONLY=ON ..
make
cd ../../

cd test_cpp
mkdir buildMakeStatic
cd buildMakeStatic
cmake -DUSE_HEADERONLY=OFF ..
make
./runUnitTests
cd ../../

