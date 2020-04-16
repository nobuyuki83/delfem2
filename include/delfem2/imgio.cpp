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
int delfem2::LoadImage_PPMAscii(
    unsigned int& width, unsigned int& height,
    std::vector<unsigned char>& image,
    const std::string& fname)
{
  FILE* fp = fopen(fname.c_str(),"r");
  if( fp == nullptr ){ return 1; }

  image.clear();
  width = height = 0;
  {
    const unsigned int buffSize = 256;
    char buff[buffSize];
    char* cres = nullptr;
    cres = fgets(buff,buffSize,fp); if( cres == nullptr ){ return 2; }
    if( buff[0] != 'P' || buff[1] != '3'){ return 2; }
    cres = fgets(buff,buffSize,fp); if( cres == nullptr ){ return 2; }
    sscanf(buff,"%d%d",&width,&height);
    cres = fgets(buff,buffSize,fp);  // read 255
    if( cres == nullptr ){ return -1; }
  }
  image.resize(width*height*3);
  const auto buffSize = (unsigned int)(4*3*width*1.2);  // ÇøÇÂÇ¡Ç∆ó]ï™ñ⁄Ç…Ç∆Ç¡ÇƒÇ®Ç≠
  std::vector<char> buff(buffSize);
  unsigned int icnt = 0;
  while (icnt<width*height*3) {
    char* cres = fgets(buff.data(),buffSize,fp);
    if( cres == nullptr ){ return -1; }
    char* pCur = buff.data();
    char* pNxt = nullptr;
    for(;;){
      //      if(      pCur[0] == ' ' ){ assert(0); }
      if(      pCur[1] == ' ' || pCur[1] == '\n' ){ pCur[1]='\0'; pNxt=pCur+2; }
      else if( pCur[2] == ' ' || pCur[2] == '\n' ){ pCur[2]='\0'; pNxt=pCur+3; }
      else if( pCur[3] == ' ' || pCur[3] == '\n' ){ pCur[3]='\0'; pNxt=pCur+4; }
      //      else{ assert(0); }
      unsigned int val = atoi(pCur);
      unsigned int ih = icnt/(width*3);
      unsigned int iw = icnt-ih*width*3;
      image[(height-ih-1)*width*3+iw] = val;
      icnt++;
      if( pNxt[0] == '\n' || pNxt[0] == '\0') break;
      pCur = pNxt;
    }
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
