echo "################################"
echo "build examples_glut"
echo "################################"

cd examples_oldgl_glut || exit
mkdir buildMake 
cd buildMake || exit
cmake ..
cmake --build .
cd ../../

cd examples_oldgl_glut || exit
mkdir buildXcode 
cd buildXcode || exit
cmake -G Xcode ..
cmake --build .
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
echo "build examples_glfwold"
echo "################################"

cd examples_oldgl_glfw || exit
mkdir buildXcodeHdronly 
cd buildXcodeHdronly || exit
cmake .. -G Xcode -DUSE_STATIC_LIB=OFF
cd ../../

cd examples_oldgl_glfw || exit
mkdir buildXcodeStatic 
cd buildXcodeStatic || exit
cmake .. -G Xcode -DUSE_STATIC_LIB=ON
cd ../../

cd examples_oldgl_glfw || exit
mkdir buildMakeHdronly
cd buildMakeHdronly || exit
cmake .. -DUSE_STATIC_LIB=OFF
cmake --build .
cd ../../

cd examples_oldgl_glfw || exit
mkdir buildMakeStatic 
cd buildMakeStatic || exit
cmake .. -DUSE_STATIC_LIB=ON
cmake --build .
cd ../../


echo "################################"
echo "build examples_glfwold_glad"
echo "################################"

cd examples_oldgl_glfw_glad || exit
mkdir buildXcode
cd buildXcode || exit
cmake .. -G Xcode 
cd ../../

cd examples_oldgl_glfw_glad || exit
mkdir buildMake
cd buildMake || exit
cmake .. 
cmake --build .
cd ../../

echo "################################"
echo "build examples_glfwnew"
echo "################################"

cd examples_newgl_glfw || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode -DUSE_STATIC_LIB=OFF
cmake --build .
cd ../../

cd examples_newgl_glfw || exit
mkdir buildXcodeStatic
cd buildXcodeStatic || exit
cmake .. -G Xcode -DUSE_STATIC_LIB=ON
cmake --build .
cd ../../

echo "################################"
echo "build examples_glfw_thread_oldgl"
echo "################################"

cd examples_oldgl_glfw_thread || exit
mkdir buildXcodeHdronly 
cd buildXcodeHdronly || exit
cmake .. -G Xcode
cmake --build .
cd ../../


echo "################################"
echo "fetch latest imgui"
echo "################################"

git submodule update --init -- 3rd_party/imgui
cd 3rd_party/imgui || exit
git checkout master
git pull origin master
cd ../../

echo "################################"
echo "build examples_newgl_glfw_imgui"
echo "################################"

cd examples_newgl_glfw_imgui || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode
cmake --build .
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
echo "compile demos using tinygltf"
echo "################################"

cd examples_oldgl_glfw_tinygltf || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode
cmake --build .
cd ../../

cd examples_oldgl_glfw_tinygltf || exit
mkdir buildXcodeStatic
cd buildXcodeStatic || exit
cmake .. -G Xcode
cmake --build .
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
echo "build examples_cnpy"
echo "################################"

cd examples_oldgl_glfw_cnpy || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode
cd ../../

cd examples_oldgl_glfw_cnpy || exit
mkdir buildMakeHdronly
cd buildMakeHdronly || exit
cmake ..
cmake --build .
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

echo "################################"
echo "build examples_eigen"
echo "################################"

cd examples_oldgl_glfw_eigen || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode
cd ../../

cd examples_oldgl_glfw_eigen || exit
mkdir buildMakeHdronly
cd buildMakeHdronly || exit
cmake ..
cmake --build .
cd ../../


echo "################################"
echo "build examples_thread"
echo "################################"

cd examples_oldgl_glfw_thread || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode
cmake --build .
cd ../../

cd examples_oldgl_glfw_thread || exit
mkdir buildMakeHdronly
cd buildMakeHdronly || exit
cmake ..
cmake --build .
cd ../../


echo "######################################"
echo "fetch latest googletest and compile it"
echo "######################################"

git submodule update --init -- 3rd_party/googletest
cd 3rd_party/googletest || exit
git checkout master
git pull origin master
cmake .
cmake --build . 
mkdir ../libgtest
cmake --install . --prefix ../libgtest
cd ../../

echo "###############################"
echo "test cpp"
echo "###############################"

cd test_cpp || exit
mkdir buildXcodeStatic
cd buildXcodeStatic || exit
cmake .. -G Xcode -DUSE_STATIC_LIB=ON
cd ../../

cd test_cpp || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode -DUSE_STATIC_LIB=OFF
cd ../../

cd test_cpp || eixt
mkdir buildMakeHdronly
cd buildMakeHdronly || eixt
cmake .. -DUSE_STATIC_LIB=OFF
cmake --build .
cd ../../

cd test_cpp || exit
mkdir buildMakeStatic
cd buildMakeStatic || eixt
cmake .. -DUSE_STATIC_LIB=ON
cmake --build .
./runUnitTests
cd ../../

echo "###############################"
echo "test eigen"
echo "###############################"

cd test_eigen || exit
mkdir buildMake
cd buildMake || exit
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
./runUnitTests
cd ../../


echo "################################"
echo "build Imath"
echo "################################"

git submodule update --init -- 3rd_party/Imath
cd 3rd_party/Imath || exit
git checkout master
git pull origin master
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
cmake .. -DALEMBIC_SHARED_LIBS=OFF -DBUILD_TESTING=OFF -DCMAKE_MODULE_PATH=../../Imathlib
cmake --build .
cmake --install . --prefix ../alembiclib
cd ../..

echo "################################"
echo "build alembic example"
echo "################################"

cd examples_alembic || exit
mkdir buildMake 
cd buildMake || exit
cmake ..
cmake --build .
cd ../../

cd examples_alembic || exit
mkdir buildXcode 
cd buildXcode || exit
cmake -G Xcode ..
cmake --build .
cd ../../
