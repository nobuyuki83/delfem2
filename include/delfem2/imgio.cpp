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

#include "delfem2/imgio.h"

// https://en.wikipedia.org/wiki/Netpbm_format
bool delfem2::LoadImage_PPMBinary
(const std::string& filename,
 std::vector<unsigned char>& image,
 int& width, int& height)
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
    image[(height-ih-1)*width*3+iw*3+ich] = ival;
    icnt++;
  }
  return 0;
}


bool delfem2::LoadTGAFile(
    const std::string& filename,
    SFile_TGA *tgaFile)
{
  FILE *filePtr;
  unsigned char ucharBad;
  short int sintBad;
  long imageSize;
  int colorMode;
  unsigned char colorSwap;
  
  // Open the TGA file.
  filePtr = fopen(filename.c_str(), "rb");
  if (filePtr == nullptr)
  {
    return false;
  }
  
  // Read the two first bytes we don't need.
  size_t n0;
  n0 = fread(&ucharBad, sizeof(unsigned char), 1, filePtr); if( n0 != 1 ){ return false; }
  n0 = fread(&ucharBad, sizeof(unsigned char), 1, filePtr); if( n0 != 1 ){ return false; }
  
  // Which type of image gets stored in imageTypeCode.
  n0 = fread(&tgaFile->imageTypeCode, sizeof(unsigned char), 1, filePtr);
  if( n0 != 1 ){ return false; }
  
  // For our purposes, the type code should be 2 (uncompressed RGB image)
  // or 3 (uncompressed black-and-white images).
  if (tgaFile->imageTypeCode != 2 && tgaFile->imageTypeCode != 3){
    fclose(filePtr);
    return false;
  }
  
  // Read 13 bytes of data we don't need.
  n0 = fread(&sintBad, sizeof(short int), 1, filePtr);  if( n0 != 1 ){ return false; }
  n0 = fread(&sintBad, sizeof(short int), 1, filePtr);  if( n0 != 1 ){ return false; }
  n0 = fread(&ucharBad, sizeof(unsigned char), 1, filePtr);  if( n0 != 1 ){ return false; }
  n0 = fread(&sintBad, sizeof(short int), 1, filePtr);  if( n0 != 1 ){ return false; }
  n0 = fread(&sintBad, sizeof(short int), 1, filePtr);  if( n0 != 1 ){ return false; }
  
  // Read the image's width and height.
  n0 = fread(&tgaFile->imageWidth, sizeof(short int), 1, filePtr);  if( n0 != 1 ){ return false; }
  n0 = fread(&tgaFile->imageHeight, sizeof(short int), 1, filePtr);  if( n0 != 1 ){ return false; }
  
  // Read the bit depth.
  n0 = fread(&tgaFile->bitCount, sizeof(unsigned char), 1, filePtr);  if( n0 != 1 ){ return false; }
  
  // Read one byte of data we don't need.
  n0 = fread(&ucharBad, sizeof(unsigned char), 1, filePtr);   if( n0 != 1 ){ return false; }
  
  // Color mode -> 3 = BGR, 4 = BGRA.
  colorMode = tgaFile->bitCount / 8;
  imageSize = tgaFile->imageWidth * tgaFile->imageHeight * colorMode;
  
  // Allocate memory for the image data.
  tgaFile->imageData = (unsigned char*)malloc(sizeof(unsigned char)*imageSize);
  
  // Read the image data.
  n0 = fread(tgaFile->imageData, sizeof(unsigned char), imageSize, filePtr);
  if( (long)n0 != imageSize ){ return false; }
  
  // Change from BGR to RGB so OpenGL can read the image data.
  for (int imageIdx = 0; imageIdx < imageSize; imageIdx += colorMode)
  {
    colorSwap = tgaFile->imageData[imageIdx];
    tgaFile->imageData[imageIdx] = tgaFile->imageData[imageIdx + 2];
    tgaFile->imageData[imageIdx + 2] = colorSwap;
  }
  
  fclose(filePtr);
  return true;
}

void delfem2::ImageInterpolation_Bilinear(
    std::vector<double>& aColor,
    int width,
    int height,
    const unsigned char* img,
    const double* aXY,
    unsigned int nXY)
{
 aColor.resize(nXY*3);
  for(unsigned int ip=0;ip<nXY;++ip){
    double x = aXY[ip*2+0]*(width-1);
    double y = (1.0-aXY[ip*2+1])*(height-1);
    int ix0 = static_cast<int>(floor(x));
    int iy0 = static_cast<int>(floor(y));
    int ix1 = ix0+1; if( ix1 == width ){ ix1 = width-1; }
    int iy1 = iy0+1; if( iy1 == height ){ iy1 = height-1; }
    assert( ix0 >= 0 && ix0 < width );
    assert( iy0 >= 0 && iy0 < height );
    assert( ix1 >= 0 && ix1 < width );
    assert( ix1 >= 0 && ix1 < width );
    double rx = x-ix0;
    double ry = y-iy0;

    double w00 = 1.0/255.0*(1-rx)*(1-ry);
    double w01 = 1.0/255.0*(1-rx)*ry;
    double w10 = 1.0/255.0*rx*(1-ry);
    double w11 = 1.0/255.0*rx*ry;
    for(int i=0;i<3;++i) {
      aColor[ip*3+i] =
          + w00 * img[(ix0 + iy0 * width) * 3 + i]
          + w01 * img[(ix0 + iy1 * width) * 3 + i]
          + w10 * img[(ix1 + iy0 * width) * 3 + i]
          + w11 * img[(ix1 + iy1 * width) * 3 + i];
    }
  }
}