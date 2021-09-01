/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_FILE_H
#define DFM2_FILE_H

#include <fstream>
#include <map>
#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

// -----------------
// Path related

DFM2_INLINE bool isFileExists(const std::string& fpath);
DFM2_INLINE std::string pathRemoveExtension(const std::string& fpath);
DFM2_INLINE std::string pathGetExtension(const std::string& fpath);
DFM2_INLINE std::string getPathDir(const std::string& fpath);


// ---------------
// File related

DFM2_INLINE bool ReadParam(
    std::vector<float>& aPara,
    const std::string& fname);


DFM2_INLINE bool GetFileContents(
    std::vector<char>& aC,
    const std::string& fpath);

template <typename T>
DFM2_INLINE bool WriteParam(
    const std::string& fname,
    const std::vector<T>& aPara)
{
  std::ofstream fout;
  fout.open(fname.c_str());
  if( !fout.is_open() ) return false;
  for(unsigned int ip=0;ip<aPara.size();++ip){
    fout << aPara[ip] << std::endl;
  }
  return true;
}

DFM2_INLINE std::string LoadFile(
    const std::string& fname);

//DFM2_INLINE std::map<std::string, std::string> ReadDictionary(
//    const std::string& path);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/file.cpp"
#endif

#endif
