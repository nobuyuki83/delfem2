
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

cd docs/00_openwin
em++ ../../examples_newgl_glfw/00_openwin/main.cpp -o index.html -s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1 -std=c++1z
cd ../../

cd docs/01_drawrect
em++  ../../examples_newgl_glfw/01_drawrect/main.cpp -o index.html -std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include"  -s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1 
cd ../../

cd docs/02_nav3d
em++  ../../examples_newgl_glfw/02_nav3d/main.cpp -o index.html -std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include"  -s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1 
cd ../../
