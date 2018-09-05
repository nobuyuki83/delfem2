#!bin/bash

echo $(uname)

if test $(uname) = "Linux"
then
  OS='LINUX'
else
  OS='OSX'
fi

echo $OS

cd external
if ! test -e "fbx"; then
  mkdir fbx
fi
cd fbx

if test $OS = 'LINUX'; then
  wget http://download.autodesk.com/us/fbx/2019/2019.0/fbx20190_fbxsdk_linux.tar.gz
  tar -xvzf fbx20190_fbxsdk_linux.tar.gz
  chmod ugo+x fbx20190_fbxsdk_linux
  ./fbx20190_fbxsdk_linux .  
elif test $OS = 'OSX'; then
  brew install wget
  wget http://download.autodesk.com/us/fbx/2019/2019.0/fbx20190_fbxsdk_clang_mac.pkg.tgz
fi


