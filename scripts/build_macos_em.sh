
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

echo "00_openwin"
cd docs/00_openwin
em++ ../../examples_newgl_glfw/00_openwin/main.cpp -o index.html \
	-std=c++1z \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1 
cd ../../

echo "01_drawrect"
cd docs/01_drawrect
em++  ../../examples_newgl_glfw/01_drawrect/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include"  \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1 
cd ../../

echo "02_nav3d"
cd docs/02_nav3d
em++  ../../examples_newgl_glfw/02_nav3d/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include" \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1 
cd ../../

echo "03_texture"
cd docs/03_texture
em++  ../../examples_newgl_glfw/03_texture/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include"  \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1 
cd ../../

echo "04_render2texture"
cd docs/04_render2texture
em++  ../../examples_newgl_glfw/04_render2texture/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include" \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1 
cd ../../

echo "05_offscreenprojection"
cd docs/05_offscreenprojection
em++  ../../examples_newgl_glfw/05_OffScreenProjection/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include" \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1
cd ../../

echo "20_cad2d"
cd docs/20_cad2d
em++  ../../examples_newgl_glfw/20_cad2d/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include" \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1
cd ../../

echo "40_femcloth"
cd docs/40_femcloth
em++  ../../examples_newgl_glfw/40_femcloth/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include" \
	-O2 -s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1
cd ../../

echo "41_fem2d_poisson"
cd docs/41_fem2d_poisson
em++  ../../examples_newgl_glfw/41_fem2d_poisson/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include" \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1
cd ../../

echo "42_fem2d_linearsolid"
cd docs/42_fem2d_linearsolid
em++  ../../examples_newgl_glfw/42_fem2d_linearsolid/main.cpp -o index.html \
	-std=c++11 -DDFM2_HEADER_ONLY=ON  -I"../../include" \
	-s USE_WEBGL2=1 -s USE_GLFW=3 -s WASM=1
cd ../../

