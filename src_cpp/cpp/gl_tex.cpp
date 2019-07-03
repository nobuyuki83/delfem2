/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdlib>

#include "delfem2/gl_tex.h"


#if defined(__APPLE__) && defined(__MACH__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else // linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif

void CTexture::LoadTex()
{
  if( id_tex == 0 ){
    ::glGenTextures(1, &id_tex);
  }
  glBindTexture(GL_TEXTURE_2D, id_tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  assert( (int)aRGB.size() == w*h*3 );
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
               w, h, 0, GL_RGB, GL_UNSIGNED_BYTE,
               aRGB.data() );
  glBindTexture(GL_TEXTURE_2D, 0);
}


void CTexture::Draw(){
  if( id_tex == 0 ){ return; }
  /*
   const CVector3& dx = x_axis;
   const CVector3& dy = Cross(z_axis,dx);
   const double lx = lengrid*nResX;
   const double ly = lengrid*nResY;
   CVector3 p0 = origin;
   CVector3 p1 = origin + lx*dx;
   CVector3 p2 = origin + lx*dx + ly*dy;
   CVector3 p3 = origin + ly*dy;
   */
  ::glEnable(GL_TEXTURE_2D);
  ::glDisable(GL_LIGHTING);
  ::glBindTexture(GL_TEXTURE_2D, id_tex);
  ::glColor3d(1,1,1);
  ::glBegin(GL_QUADS);
  ::glTexCoord2d(0.0, 0.0); ::glVertex3d(0,0,0);
  ::glTexCoord2d(1.0, 0.0); ::glVertex3d(w,0,0);
  ::glTexCoord2d(1.0, 1.0); ::glVertex3d(w,h,0);
  ::glTexCoord2d(0.0, 1.0); ::glVertex3d(0,h,0);
  ::glEnd();
  ::glBindTexture(GL_TEXTURE_2D, 0);
  ::glDisable(GL_TEXTURE_2D);
}


/////////////////////////////////////////////////////////////////////////////////////

void SaveImage(const std::string& path)
{
  static unsigned int inum = 0;
  int viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  void* image = malloc(3*viewport[2]*viewport[3]);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, viewport[2], viewport[3], GL_RGB, GL_UNSIGNED_BYTE, image);
  unsigned int width = viewport[2];
  unsigned int height = viewport[3];
  std::ofstream fout;
  //  //  fname << "out";
  fout.open(path.c_str(), std::ios::out);
  fout<<"P3\n";
  fout<<width<<" "<<height<<"\n";
  fout<<"255\n";
  //  fout << "255\n";
  //  fout << "255\n";
  char* img = (char*)image;
  for (unsigned int ih = 0; ih<height; ih++){
    for (unsigned int iw = 0; iw<width; iw++){
      unsigned int i = (height-1-ih)*width+iw;
      int r = (unsigned char)img[i*3+0];
      int g = (unsigned char)img[i*3+1];
      int b = (unsigned char)img[i*3+2];
      fout<<r<<" "<<g<<" "<<b<<"\n";
      //    std::cout << i << " " << r << " "<< g << " "<< b << std::endl;
    }
  }
  fout.close();
  //  if( inum >= 700 ) abort();
  //  if( inum >= 400 ) abort();
  if (inum>=600) abort();
  inum++;
}


void LoadImage_PPM
(const std::string& filename,
 std::vector<unsigned char>& image,
 int& width, int& height)
{
  std::ifstream file(filename.c_str(), std::ios::binary);
  if (!file) {
    std::cerr<<"Could not open file \""<<filename<<"\"."<<std::endl;
    return;
  }
  
  std::string header;
  {
    char buff[256];
    file.getline(buff, 256);
    header = std::string(buff);
  }
  if (header=="P6") {
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
}

int ReadPPM_SetTexture(const std::string& fname)
{
  std::cout << "ReadPPM " << std::endl;
  FILE* fp = fopen(fname.c_str(),"r");
  if( fp == NULL ){
    std::cout << "Read PPM Fail" << std::endl;
    return -1;
  }
  
  int w, h;
  std::vector<char> aRGB;
  {
    const unsigned int buffSize = 256;
    char buff[buffSize];
    fgets(buff,buffSize,fp); std::cout << buff << std::endl;
    fgets(buff,buffSize,fp); std::cout << buff << std::endl;
    sscanf(buff,"%d%d",&w,&h);
    fgets(buff,buffSize,fp);  // read 255
  }
  std::cout << "tex size : " << w << " " << h << std::endl;
  //  assert( w >= 0 && h >=0 );
  aRGB.resize(w*h*3);
  const unsigned int buffSize = (unsigned int)(4*3*w*1.2);  // ÇøÇÂÇ¡Ç∆ó]ï™ñ⁄Ç…Ç∆Ç¡ÇƒÇ®Ç≠
  char* buff = new char [buffSize];
  int icnt = 0;
  while (icnt<w*h*3) {
    fgets(buff,buffSize,fp);
    char* pCur = buff;
    char* pNxt;
    for(;;){
      //      if(      pCur[0] == ' ' ){ assert(0); }
      if(      pCur[1] == ' ' || pCur[1] == '\n' ){ pCur[1]='\0'; pNxt=pCur+2; }
      else if( pCur[2] == ' ' || pCur[2] == '\n' ){ pCur[2]='\0'; pNxt=pCur+3; }
      else if( pCur[3] == ' ' || pCur[3] == '\n' ){ pCur[3]='\0'; pNxt=pCur+4; }
      //      else{ assert(0); }
      unsigned int val = atoi(pCur);
      unsigned int ih = icnt/(w*3);
      unsigned int iw = icnt-ih*w*3;
      aRGB[(h-ih-1)*w*3+iw] = val;
      icnt++;
      if( pNxt[0] == '\n' || pNxt[0] == '\0') break;
      pCur = pNxt;
    }
  }
  delete[] buff;
  //  this->SetImage(w,h,aRGB);
  
  std::cout << "width height : " << w << " " << h << std::endl;
  
  ////////////////
  
  GLubyte* inputRGB = new GLubyte [w*h*3];
  for(int i=0;i<w*h*3;i++){ inputRGB[i] = aRGB[i]; }
  
  int m_texWidth = 256;
  int m_texHeight = 256;
  //  std::cout << m_texWidth << " " << m_texHight << std::endl;
  
  GLubyte* scaledRGB;
  if( w == m_texWidth && h == m_texHeight ){
    scaledRGB = inputRGB;
  }
  else{
    scaledRGB = new GLubyte [m_texWidth*m_texHeight*3];
    gluScaleImage( GL_RGB, w, h, GL_UNSIGNED_BYTE, inputRGB,
                  m_texWidth, m_texHeight, GL_UNSIGNED_BYTE, scaledRGB );
    delete [] inputRGB;
  }
  
  glEnable(GL_TEXTURE_2D);
  GLuint m_texName = 0;
  glGenTextures(1 , &m_texName);
  glBindTexture(GL_TEXTURE_2D , m_texName);
  glTexImage2D(GL_TEXTURE_2D , 0 , 3 , m_texWidth, m_texHeight,
               0 , GL_RGB , GL_UNSIGNED_BYTE , scaledRGB );
  delete[] scaledRGB;
  
  //  std::cout << m_texName << std::endl;
  
  return (int)m_texName;
}

bool LoadTGAFile
(const char *filename,
 SFile_TGA *tgaFile)
{
  FILE *filePtr;
  unsigned char ucharBad;
  short int sintBad;
  long imageSize;
  int colorMode;
  unsigned char colorSwap;
  
  // Open the TGA file.
  filePtr = fopen(filename, "rb");
  if (filePtr == NULL)
  {
    return false;
  }
  
  // Read the two first bytes we don't need.
  fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
  fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
  
  // Which type of image gets stored in imageTypeCode.
  fread(&tgaFile->imageTypeCode, sizeof(unsigned char), 1, filePtr);
  
  // For our purposes, the type code should be 2 (uncompressed RGB image)
  // or 3 (uncompressed black-and-white images).
  if (tgaFile->imageTypeCode != 2 && tgaFile->imageTypeCode != 3)
  {
    fclose(filePtr);
    return false;
  }
  
  // Read 13 bytes of data we don't need.
  fread(&sintBad, sizeof(short int), 1, filePtr);
  fread(&sintBad, sizeof(short int), 1, filePtr);
  fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
  fread(&sintBad, sizeof(short int), 1, filePtr);
  fread(&sintBad, sizeof(short int), 1, filePtr);
  
  // Read the image's width and height.
  fread(&tgaFile->imageWidth, sizeof(short int), 1, filePtr);
  fread(&tgaFile->imageHeight, sizeof(short int), 1, filePtr);
  
  // Read the bit depth.
  fread(&tgaFile->bitCount, sizeof(unsigned char), 1, filePtr);
  
  // Read one byte of data we don't need.
  fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
  
  // Color mode -> 3 = BGR, 4 = BGRA.
  colorMode = tgaFile->bitCount / 8;
  imageSize = tgaFile->imageWidth * tgaFile->imageHeight * colorMode;
  
  // Allocate memory for the image data.
  tgaFile->imageData = (unsigned char*)malloc(sizeof(unsigned char)*imageSize);
  
  // Read the image data.
  fread(tgaFile->imageData, sizeof(unsigned char), imageSize, filePtr);
  
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

GLuint LoadTexture
(const unsigned char* image,
 const int width, const int height, const int bpp)
{
  GLuint id_tex; glGenTextures(1, &id_tex);
  glBindTexture(GL_TEXTURE_2D, id_tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  if( bpp == 3 ){
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, (GLsizei)width, (GLsizei)height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
  }
  if( bpp == 4 ){
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, (GLsizei)width, (GLsizei)height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
  }
  return id_tex;
};

void DrawTextureBackground
(const GLuint tex,
 const int imgWidth,
 const int imgHeight,
 const int winWidth,
 const int winHeight)
{
  double imgAsp = (double)imgWidth/imgHeight;
  double winAsp = (double)winWidth/winHeight;
  /////
  glPushAttrib(GL_TRANSFORM_BIT|GL_CURRENT_BIT|GL_ENABLE_BIT);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  
  if (winAsp>imgAsp){
    double imgWidth2 = imgHeight*winAsp;
    gluOrtho2D(
               -0.5*imgWidth2,+0.5*imgWidth2,
               -0.5*imgHeight,+0.5*imgHeight);
  }
  else{
    double imgHeight2 = (double)imgWidth/winAsp;
    gluOrtho2D(
               -0.5*imgWidth, +0.5*imgWidth,
               -0.5*imgHeight2, +0.5*imgHeight2);
  }
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  
  glBindTexture(GL_TEXTURE_2D, tex);
  
  glColor3f(1.0f, 1.0f, 1.0f);
  glBegin(GL_QUADS);
  glTexCoord2i(0, 0); glVertex2d(-0.5*imgWidth, -0.5*imgHeight);
  glTexCoord2i(1, 0); glVertex2d(+0.5*imgWidth, -0.5*imgHeight);
  glTexCoord2i(1, 1); glVertex2d(+0.5*imgWidth, +0.5*imgHeight);
  glTexCoord2i(0, 1); glVertex2d(-0.5*imgWidth, +0.5*imgHeight);
  glEnd();
  
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  
  glPopAttrib();
}

// use it for GLSL shader drawing
void DrawRectangle_FullCanvas()
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  
  glBegin(GL_QUADS);
  glVertex2d(-1, -1);
  glVertex2d(+1, -1);
  glVertex2d(+1, +1);
  glVertex2d(-1, +1);
  glEnd();
  
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  
  glPopAttrib();
}


void CTexManager::Clear(){
  for(int itex=0;itex<(int)aTexInfo.size();++itex){
    unsigned int id_tex_gl = aTexInfo[itex].id_tex_gl;
    if( glIsTexture(id_tex_gl) ){
      ::glDeleteTextures(1, &id_tex_gl);
    }
  }
  aTexInfo.clear();
}

void CTexManager::BindTexturePath(const std::string& path) const {
  for(int itex=0;itex<(int)aTexInfo.size();++itex){
    if( aTexInfo[itex].full_path != path ) continue;
    glBindTexture(GL_TEXTURE_2D, aTexInfo[itex].id_tex_gl );
    glEnable(GL_TEXTURE_2D);
  }
}
