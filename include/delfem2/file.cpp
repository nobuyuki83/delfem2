/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <sstream>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include "delfem2/file.h"

namespace delfem2{
namespace file{

DFM2_INLINE void Split(
    std::vector<std::string>& aToken,
    const std::string& str,
    char delimiter)
{
  aToken.clear();
  std::stringstream data(str);
  std::string line;
  while(std::getline(data,line,delimiter)){
    if( line.empty() ){ continue; }
    aToken.push_back(line);
  }
}

}
}


// -------------------------------------
// file IO

DFM2_INLINE std::string delfem2::LoadFile
(const std::string& fname)
{
  std::ifstream inputFile1;
  inputFile1.open(fname.c_str());
  if( !inputFile1.is_open() ){
    std::cout << "Error! --> cannot open the file: " << fname << std::endl;
    return "";
  }
  std::istreambuf_iterator<char> vdataBegin(inputFile1);
  std::istreambuf_iterator<char> vdataEnd;
  return std::string(vdataBegin,vdataEnd);
}

DFM2_INLINE bool delfem2::isFileExists(const std::string& fpath)
{
  std::ifstream fin;
  fin.open(fpath.c_str());
  return fin.is_open();
}

DFM2_INLINE std::string delfem2::pathRemoveExtension(const std::string& fpath)
{
  std::vector<std::string> aToken;
  file::Split(aToken, fpath, '.');
  std::string sRes;
  for(int it=0;it<(int)aToken.size()-1;++it){
    sRes += aToken[it];
  }
  return sRes;
}

DFM2_INLINE std::string delfem2::pathGetExtension(const std::string& fpath)
{
  std::vector<std::string> aToken;
  file::Split(aToken, fpath, '.');
  std::string sRes;
  if( !aToken.empty() ){
    sRes = aToken[aToken.size()-1];
  }
  return sRes;
}

DFM2_INLINE std::string delfem2::getPathDir(const std::string& fpath)
{
  const size_t iloc = fpath.find_last_of('/');
  std::string sres = std::string(fpath.begin(),fpath.begin()+iloc);
  return sres;
}

// ----------------------------------

DFM2_INLINE bool delfem2::ReadParam(
    std::vector<float>& aPara,
    const std::string& fname)
{
  std::ifstream fin;
  fin.open(fname.c_str());
  if( !fin.is_open() ) return false;
  aPara.clear();
  for(;;){
    double d;
    fin >> d;
    if( fin.eof() ) break;
    aPara.push_back(static_cast<float>(d));
  }
  return true;
}

// ----------------------

DFM2_INLINE bool delfem2::GetFileContents(
    std::vector<char>& aC,
    const std::string& fpath)
{
  FILE* fp = nullptr;
  size_t size;
  
  fp = fopen(fpath.c_str(), "rb");
  if (!fp) goto error;
  fseek(fp, 0, SEEK_END);
  size = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  aC.resize(size+1);
  if (fread(aC.data(), 1, size, fp) != size) goto error;
  aC[size] = '\0';  // Must be null terminated.
  fclose(fp);
  return true;
  
error:
  if (fp) fclose(fp);
  return false;
}