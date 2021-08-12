/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/filenpy_str.h"

#include <map>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include "delfem2/str.h"

namespace delfem2{
namespace funcs{

// 10bytes header
struct NPY
{
  char magic_string[6] = {'X','X','X','X','X','X'}; // 6 bytes (0x93NUMPY)
  unsigned char major_version = 0; // 1 byte
  unsigned char minor_version = 0; // 1 byte
  unsigned short header_len = 0; // 2 bytes
};

DFM2_INLINE bool LoadNumpy(
    int& ndim0,
    int& ndim1,
    std::ifstream& fp)
{
  NPY npy;
//  size_t n0 = fread(&npy, sizeof(npy), 1, fp);
//  if( n0 != 1 ){ return false; }
  fp.read((char*)&npy, sizeof(npy));
  if( fp.fail() ){ return false; }

  { // check magic string
    unsigned char sMagic[7] = {0x93,'N','U','M','P','Y'};
    if( memcmp(npy.magic_string, sMagic, 6 ) != 0 ){ return false; }
  }

  ndim0 = 0;
  ndim1 = 0;
  { // format
    char buff[256];
//    size_t n1 = fread(buff, 1, npy.header_len, fp);
//    if( n1 != npy.header_len ){ return false; }
    fp.read(buff,npy.header_len);
    if( fp.fail() ){ return false; }
    std::map<std::string, std::string> map0 = ReadDictionary_Python(std::string(buff));
    std::string str_shape = map0["'shape'"];
    str_shape = Get_Parentheses(str_shape,"()");
    std::vector<std::string> aToken = Split(str_shape,',');
    if( aToken.size() != 2 ){ return false; }
    ndim0 = myStoi(aToken[0]);
    ndim1 = myStoi(aToken[1]);
  }
  return true;
}


}
}


// ------------------------------
// -----------------------------

// Read somehting like this {'descr': '<f8', 'fortran_order': False, 'shape': (3682, 151), }
std::map<std::string, std::string>
DFM2_INLINE delfem2::ReadDictionary_Python(
    const std::string& strIn)
{
  std::string buff = RemoveSpace(strIn);
  std::map<std::string, std::string> map0;
  //
  const size_t n =  buff.length();
  const char* pbuff = buff.c_str();
  //
  int iend = 1; // first letter is {
  int icolon = -1;
  bool is_parentheses = false;
  std::string skey;
  for(unsigned int i=0;i<n;++i){
    if( buff[i] == ':' ){
      char str[256];
//      strncpy(str, pbuff+iend,i-iend);
      for(unsigned int j=0;j<i-iend;++j){ str[j] = pbuff[iend+j]; }
      str[i-iend] = '\0';
      //      std::cout << "key:" << str << "#" << std::endl;
      skey = std::string(str);
      //      if( buff[i+1] == ' ' ){ i+=1; }
      icolon = i;
    }
    if( buff[i] == ','){
      assert( icolon != -1 );
      if( is_parentheses ){ continue; }
      char str[256];
//      strncpy(str, pbuff+icolon+1,i-icolon-1);
      for(unsigned int j=0;j<i-icolon-1;++j){ str[j] = pbuff[icolon+1+j]; }
      str[i-icolon-1] = '\0';
      //      std::cout << "val:" << str << "#" << std::endl;
      std::string svalue(str);
      map0.insert( std::make_pair(skey,svalue) );
      //      if( buff[i+1] == ' ' ){ i+=1; }
      iend = i+1;
      icolon = -1;
    }
    if( buff[i] == '('){
      is_parentheses = true;
    }
    if( buff[i] == ')'){
      is_parentheses = false;
    }
  }
  return map0;
}

//bool isNaN(double x) { return x!=x; }

template <typename REAL>
DFM2_INLINE bool delfem2::LoadNumpy_2Dim(
    int& ndim0,
    int& ndim1,
    std::vector<REAL>& aData,
    const std::string& path)
{
//  FILE* fp = fopen(path.c_str(),"rb");
  std::ifstream fin(path,std::ios::in | std::ios::binary);
  if( fin.fail() ) { return false; }
  funcs::LoadNumpy(ndim0, ndim1, fin);
  int size = ndim0*ndim1;
  aData.resize(size);
//  size_t n0 = fread(&aData[0], sizeof(float), size, fp);
//  return (int) n0 == size;
  fin.read( (char*)&aData[0], sizeof(REAL)*size );
  return !fin.fail();
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::LoadNumpy_2Dim(
    int&, int&, std::vector<float>& aData, const std::string& path);
template bool delfem2::LoadNumpy_2Dim(
    int&, int&, std::vector<double>& aData, const std::string& path);
#endif



DFM2_INLINE bool delfem2::LoadNumpy_1DimF(
    int& ndim0,
    std::vector<float>& aData,
    const std::string& path)
{
//  FILE* fp = fopen(path.c_str(),"r");
//  if( fp == nullptr ) { return false; }

  std::ifstream fin(path,std::ios::in | std::ios::binary);
  if( fin.fail() ) { return false; }
  
  funcs::NPY npy;
//  size_t n0 = fread(&npy, sizeof(npy), 1, fp);
//  if( n0 != 1 ){ return false; }
  fin.read((char*)&npy,sizeof(npy));
  if( fin.fail() ){ return false; }
  
  { // check magic string
    unsigned char sMagic[7] = {0x93,'N','U','M','P','Y'};
    if( memcmp(npy.magic_string, sMagic, 6 ) != 0 ){ return false; }
  }
  
  //    std::cout << npy.major_version << std::endl;
  //    std::cout << npy.minor_version << std::endl;
  //    std::cout << npy.header_len << std::endl;
  
  ndim0 = 0;
  { // format
    char buff[256];
//    size_t n1 = fread(buff, 1, npy.header_len, fp);
//    if( n1 != npy.header_len ){ return false; }
    fin.read(buff,npy.header_len);
    if(fin.fail()){ return false; }
    //    std::cout << buff << "###" << strlen(buff) << std::endl;
    std::map<std::string, std::string> map0 = ReadDictionary_Python(std::string(buff));
    //    for(std::map<std::string, std::string>::iterator itr=map0.begin();itr!=map0.end();++itr){
    //      std::cout << itr->first << " --> " << itr->second << std::endl;
    //    }
    std::string str_shape = map0["'shape'"];
    str_shape = Get_Parentheses(str_shape,"()");
    std::vector<std::string> aToken = Split(str_shape,',');
    if( aToken.size() != 1 ){ return false; }
    ndim0 = myStoi(aToken[0]);
  }
  //    std::cout << ndim0 << " " << ndim1 << std::endl;
  
  //////
  unsigned int size = ndim0;
  aData.resize(size);
  //  double* aRes = (double*)malloc( sizeof(double)*size );
//  size_t n2 = fread(&aData[0], sizeof(float), size, fp);
//  return n2 == size;
  fin.read( (char*)&aData[0], sizeof(float)*size );
  return !fin.fail();
}
