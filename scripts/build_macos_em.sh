
echo "################################"
echo "fetch & compile glfw"
echo "################################"

git submodule update --init -- 3rd_party/glfw
cd 3rd_party/glfw
git checkout master
git pull origin master
cmake .
make
cd ../..

echo "################################"
echo "fetch imgui"
echo "################################"

git submodule update --init -- 3rd_party/imgui
cd 3rd_party/imgui
git checkout master
git pull origin master
cd ../../

echo "################################"
echo "build examples_glfwnew"
echo "################################"

cd examples_glfwnew
mkdir buildEm
cd buildEm
cmake -DEMSCRIPTEN=ON -DUSE_HEADERONLY=ON ..
make
cd ../../

