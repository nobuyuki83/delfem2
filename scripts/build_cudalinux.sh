
echo "#############################"
echo "## download&build googletest"

git submodule update --init -- 3rd_party/googletest
cd 3rd_party/googletest
git checkout master 
git pull origin master
cmake .
cmake --build .
cmake --install . --prefix ../libgtest
cd ../..

echo "#######################"
echo "## cuda                "

cd examples_cuda
mkdir buildMake
cd buildMake
cmake ..
make
cd ../../

cd test_cuda
mkdir buildMake
cd buildMake
cmake ..
make
./test_cuda
cd ../../

