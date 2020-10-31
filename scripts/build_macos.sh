echo "#################################"
echo "download & build submodules"
echo "#################################"

git submodule update --init --recursive
git submodule foreach git pull origin master


echo "################################"
echo "build examples_glfwnew"
echo "################################"

cd 3rd_party/glfw
cmake .
make
cd ../..

cd examples_glfwnew
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
cmake --build .
cd ../../

cd examples_glfwnew
mkdir buildXcodeStatic
cd buildXcodeStatic
cmake -G Xcode -DUSE_HEADERONLY=OFF ..
cmake --build .
cd ../../

cd examples_glfwnew
mkdir buildEm
cd buildEm
cmake -DEMSCRIPTEN=ON -DUSE_HEADERONLY=ON ..
make
cd ../../


echo "################################"
echo "build examples_glfwold"
echo "################################"

cd examples_glfwold
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
# cmake --build . # skip build to save time
cd ../../

cd examples_glfwold
mkdir buildXcodeStatic
cd buildXcodeStatic
cmake -G Xcode -DUSE_HEADERONLY=OFF ..
# cmake --build . # skip build to save time
cd ../../

cd examples_glfwold
mkdir buildMakeHdronly 
cd buildMakeHdronly
cmake -DUSE_HEADERONLY=ON ..
make
cd ../../

cd examples_glfwold
mkdir buildMakeStatic 
cd buildMakeStatic
cmake -DUSE_HEADERONLY=OFF ..
make
cd ../../


echo "################################"
echo "build examples_smpl"
echo "################################"

cd examples_smpl
mkdir buildXcodeHdronly
cd buildXcodeHdronly
cmake -G Xcode -DUSE_HEADERONLY=ON ..
# cmake --build . # skip build to save time
cd ../../

cd examples_smpl
mkdir buildMakeHdronly
cd buildMakeHdronly
cmake -DUSE_HEADERONLY=ON ..
make 
cd ../../


echo "################################"
echo "build examples_glut"
echo "################################"

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


echo "###############################"
echo "test cpp"
echo "###############################"
# (this takes time so put it in the end)

cd "3rd_party/googletest"
cmake .
make
cd ../../

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


echo "################################"
echo "# SMPL"
echo "################################"
pip3 install chumpy
cd test_inputs
python3 smpl_preprocess.py 
pip3 uninstall chumpy -y


