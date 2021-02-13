echo "#######################"
echo "# glut"

cd examples_oldgl_glut
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../


echo "#############################"
echo "## download&build googletest"

git submodule update --init -- 3rd_party/googletest
cd 3rd_party/googletest
git checkout master 
git pull origin master
cmake .
make
cd ../..

echo "#############################"
echo "## download tinygltf"

git submodule update --init -- 3rd_party/tinygltf
cd 3rd_party/tinygltf
git checkout master 
git pull origin master
cd ../..

cd test_cpp
mkdir buildMake
cd buildMake
cmake ..
make
./runUnitTests
cd ../../



echo "test for C++ finished"

#######################

git submodule update --init -- 3rd_party/glfw
cd 3rd_party/glfw
git checkout master
git pull origin master
cmake .
make
cd ../..

cd examples_oldgl_glfw
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../


git submodule update --init -- 3rd_party/imgui
cd 3rd_party/imgui
git checkout master
git pull origin master
cd ../..

cd examples_newgl_glfw
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../


