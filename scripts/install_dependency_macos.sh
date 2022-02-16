echo "################################"
echo "fetch latest stb"
echo "################################"

git submodule update --init -- 3rd_party/stb
cd 3rd_party/stb || exit
git checkout master
git pull origin master
cd ../../

echo "################################"
echo "fetch latest glfw and compile it"
echo "################################"

git submodule update --init -- 3rd_party/glfw
cd 3rd_party/glfw || exit
git checkout master
git pull origin master
cmake .
cmake --build . 
mkdir ../libglfw
cmake --install . --prefix ../libglfw
cd ../..

echo "################################"
echo "fetch latest imgui"
echo "################################"

git submodule update --init -- 3rd_party/imgui
cd 3rd_party/imgui || exit
git checkout master
git pull origin master
cd ../../

git submodule update --init -- 3rd_party/ImGuiFileDialog
cd 3rd_party/imguiFileDialog  || exit
git checkout Lib_Only
git pull origin Lib_Only
cd ../../

echo "################################"
echo "fetch latest tinygltf"
echo "################################"

git submodule update --init -- 3rd_party/tinygltf
cd 3rd_party/tinygltf || exit
git checkout master
git pull origin master
cd ../../

echo "################################"
echo "fetch latest cnpy"
echo "################################"

git submodule update --init -- 3rd_party/cnpy
cd 3rd_party/cnpy || exit
git checkout master
git pull origin master
cd ../../

echo "################################"
echo "fetch latest eigen"
echo "################################"

git submodule update --init -- 3rd_party/eigen
cd 3rd_party/eigen || exit
git checkout master
git pull origin master
mkdir build 
cd build || exit
cmake ..
cmake --install . --prefix ../../libeigen
cd ../../../

echo "######################################"
echo "fetch latest googletest and compile it"
echo "######################################"

git submodule update --init -- 3rd_party/googletest
cd 3rd_party/googletest || exit
git checkout main
git pull origin main
cmake .
cmake --build . 
mkdir ../libgtest
cmake --install . --prefix ../libgtest
cd ../../

echo "################################"
echo "build Imath"
echo "################################"

git submodule update --init -- 3rd_party/Imath
cd 3rd_party/Imath || exit
git checkout main
git pull origin main
mkdir build
cd build || exit
cmake .. -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF
cmake --build .
cmake --install . --prefix ../../Imathlib
cd ../../../

echo "################################"
echo "build alembic"
echo "################################"

git submodule update --init -- 3rd_party/alembic
cd 3rd_party/alembic || exit
git checkout master
git pull origin master
mkdir build
cd build || exit
Imath_dir=$(pwd)/../../Imathlib
echo "Imath_dir: ${Imath_dir}"
cmake .. -DALEMBIC_SHARED_LIBS=OFF -DBUILD_TESTING=OFF -DCMAKE_PREFIX_PATH=${Imath_dir}
cmake --build . --config Release
cmake --install . --prefix ../../alembiclib
cd ../../../

echo "################################"
echo "fetch latest pugixml"
echo "################################"

git submodule update --init -- 3rd_party/pugixml
cd 3rd_party/pugixml || exit
git checkout master
git pull origin master
cd ../../

