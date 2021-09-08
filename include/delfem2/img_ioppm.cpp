/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "delfem2/img_ioppm.h"

// https://en.wikipedia.org/wiki/Netpbm_format
bool delfem2::LoadImage_PPMBinary(
    const std::string& filename,
    std::vector<unsigned char>& image,
    int& width,
    int& height)
{
  std::ifstream file(filename.c_str(), std::ios::binary);
  if (!file) {
    return false;
  }
  
  std::string header;
  {
    char buff[256];
    file.getline(buff, 256);
    header = std::string(buff);
  }
  if (header=="P6") { // binary format
    {
      int max;
      char buff[256];
      file.getline(buff,256);
      if (buff[0]=='#'){
        file>>width>>height>>max;
      }
      else{
        std::stringstream ss(buff);
        ss>>width>>height>>max;
      }
      //      std::cout<<header<<" "<<width<<" "<<height<<" "<<max<<std::endl;
    }
    // max is supporse to be 255
    const int size = width*height*3+1;
    std::vector<unsigned char> data(size);
    file.read(reinterpret_cast<char*>(&(data[0])), size);
    image.resize(3*height*width+256);
    for (int row = 0; row < height; ++row) {
      for (int col = 0; col < width; ++col) {
        int dest_index = 3*((height-row-1) * width+col);
        int src_index = (row * width+col)*3;
        image[dest_index+0] = data[src_index+1];
        image[dest_index+1] = data[src_index+2];
        image[dest_index+2] = data[src_index+3];
      }
    }
  }

  file.close();
  return true;
}

// https://en.wikipedia.org/wiki/Netpbm_format
// TODO: cannot process comment out sybol yet '#'
int delfem2::LoadImage_PPMAscii(
    unsigned int& width, 
	unsigned int& height,
    std::vector<unsigned char>& image,
    const std::string& fpath)
{
  std::ifstream fin;
  fin.open(fpath.c_str(), std::ios::in);
  if( fin.fail() ){ return 1; }

  image.clear();
  width = height = 0;
  int vmax;
  {
    const unsigned int buffSize = 256;
    char buff[buffSize];
    // get format
    if( !fin.getline(buff, buffSize) ){ return 2; }
    if( buff[0] != 'P' || buff[1] != '3'){ return 2; }
    // get size
    if( !fin.getline(buff, buffSize) ){ return 2; }
    {
      std::stringstream ss(buff);
      ss >> width >> height;
    }
    if( !fin.getline(buff, buffSize) ){ return 2; }
    {
      std::stringstream ss(buff);
      ss >> vmax;
    }
  }
  image.resize(width*height*3);
  const auto buffSize = (unsigned int)(4*3*width*1.2);
  std::vector<char> buff(buffSize, 0);
  unsigned int icnt = 0;
  while (icnt<width*height*3) {
    int ival;
    fin >> ival;
    unsigned int ih = icnt/(width*3);
    unsigned int iw = (icnt-ih*width*3)/3;
    unsigned int ich = icnt - ih*width*3 - iw*3;
    image[(height-ih-1)*width*3+iw*3+ich] = static_cast<unsigned char>(ival);
    icnt++;
  }
  return 0;
}
