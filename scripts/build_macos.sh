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
cd ../../

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
mkdir buildMakeStatic
cd buildMakeStatic || exit
cmake .. -DUSE_STATIC_LIB=ON
cmake --build .
cd ../../

cd examples_newgl_glfw || exit
mkdir buildXcodeStatic
cd buildXcodeStatic || exit
cmake .. -G Xcode -DUSE_STATIC_LIB=ON
cd ../../

cd examples_newgl_glfw || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode
cd ../../

echo "################################"
echo "build examples_newgl_glfw_imgui"
echo "################################"

cd examples_newgl_glfw_imgui || exit
mkdir buildMake
cd buildMake || exit
cmake .. 
cmake --build .
cd ../../

cd examples_newgl_glfw_imgui || exit
mkdir buildXCode
cd buildXcode || exit
cmake .. -G Xcode
cd ../../


echo "################################"
echo "compile demos using tinygltf"
echo "################################"

cd examples_oldgl_glfw_tinygltf || exit
mkdir buildMake
cd buildMake || exit
cmake .. 
cmake --build .
cd ../../

cd examples_oldgl_glfw_tinygltf || exit
mkdir buildXcode
cd buildXcode || exit
cmake .. -G Xcode
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
echo "build examples_eigen"
echo "################################"

cd examples_oldgl_glfw_eigen || exit
mkdir buildXcodeHdronly
cd buildXcodeHdronly || exit
cmake .. -G Xcode
cd ../../

cd examples_oldgl_glfw_eigen || exit
mkdir buildXcodeStatic
cd buildXcodeStatic || exit
cmake .. -G Xcode -DUSE_STATIC_LIB=ON
cd ../../

cd examples_oldgl_glfw_eigen || exit
mkdir buildMakeHdronly
cd buildMakeHdronly || exit
cmake ..
cmake --build .
cd ../../

echo "###############################"
echo "test"
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
echo "test with eigen"
echo "###############################"

cd test_cpp/eigen || exit
mkdir buildMake
cd buildMake || exit
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
./runUnitTests
cd ../../../

cd test_cpp/eigen || exit
mkdir buildXcode
cd buildXcode || exit
cmake .. -G Xcode 
cmake --build . 
cd ../../../

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
cd ../../

echo "################################"
echo "build example with pugixml"
echo "################################"

cd examples_oldgl_glfw_pugixml || exit
mkdir buildMake 
cd buildMake || exit
cmake ..
cmake --build .
cd ../../

cd examples_oldgl_glfw_pugixml || exit
mkdir buildXcode
cd buildXcode || exit
cmake -G Xcode ..
cmake --build .
cd ../../
