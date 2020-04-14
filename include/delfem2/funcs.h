/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_FUNCS_H
#define DFM2_FUNCS_H

#include "delfem2/dfm2_inline.h"
#include <fstream>
#include <map>
#include <vector>
#include <string>

namespace delfem2 {

DFM2_INLINE int myStoi(const std::string& str);
DFM2_INLINE double myStod(const std::string& str);

// --------------------
// command related

DFM2_INLINE std::string getCmdOptionStr(char ** begin, char ** end, const std::string & option);
DFM2_INLINE int getCmdOptionInt(char ** begin, char ** end, const std::string & option, int ndef);
DFM2_INLINE double getCmdOptionDouble(char ** begin, char ** end, const std::string & option, int ddef);
DFM2_INLINE bool cmdOptionExists(char** begin, char** end, const std::string& option);

// -----------------
// Path related

DFM2_INLINE bool isFileExists(const std::string& fpath);
DFM2_INLINE std::string pathRemoveExtension(const std::string& fpath);
DFM2_INLINE std::string pathGetExtension(const std::string& fpath);
DFM2_INLINE std::string getPathDir(const std::string& fpath);

// ---------------
// File related

DFM2_INLINE bool ReadParam(std::vector<float>& aPara,
               const std::string& fname);

template <typename T>
DFM2_INLINE bool WriteParam
(const std::string& fname,
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

DFM2_INLINE std::string LoadFile(const std::string& fname);

// ---------------------
// string handling

DFM2_INLINE std::vector<std::string> Split(const std::string& str,
                               char delimiter);
DFM2_INLINE void Split(std::vector<std::string>& aToken,
                       const std::string& str,
                       char delimiter);
DFM2_INLINE std::vector<std::string> Split(const std::string& str,
                                           const std::string& del);
DFM2_INLINE std::vector<std::string> Split_Parentheses(const std::string& str,
                                                       char delimiter,
                                                       const std::string& par);
DFM2_INLINE std::vector<std::string> Split_Quote(const std::string& str,
                                     char delimiter,
                                     char quote);
DFM2_INLINE std::string Get_Parentheses(const std::string& str,
                                        const std::string& par);
DFM2_INLINE std::string Replace(const std::string& str,
                                const char cf, const char ct);
DFM2_INLINE std::string Remove(const std::string& str,
                               const std::string& del);
DFM2_INLINE std::string RemoveSpace(const std::string& str);
DFM2_INLINE std::string RemoveBeginning(const std::string& str,
                                        const std::string& del);

DFM2_INLINE std::map<std::string, std::string> ReadDictionary(const std::string& path);
DFM2_INLINE std::map<std::string, std::string> ReadDictionary_Json(const std::string& strIn);

DFM2_INLINE void ReadVector_CSV(std::vector<double>& str, const std::string& strIn);

// ------------------------------
// Python realted funcs

DFM2_INLINE std::map<std::string, std::string> ReadDictionary_Python(const std::string& buff);

DFM2_INLINE bool LoadNumpy_2DimF
 (int& ndim0, int& ndim1, std::vector<float>& aData,
  const std::string& path);

DFM2_INLINE bool LoadNumpy_2DimD
 (int& ndim0, int& ndim1, std::vector<double>& aData,
  const std::string& path);

DFM2_INLINE bool LoadNumpy_1DimF
 (int& ndim0, std::vector<float>& aData,
  const std::string& path);



// ---------------------
// SVG related function

DFM2_INLINE std::string Str_SVGPolygon
 (const std::vector<double>& aXY,
  double scale);


// ---------------------------------------------------

DFM2_INLINE bool GetFileContents
 (std::vector<char>& aC,
  const std::string& fpath);

/**
 * @brief make associative array
 * @param[in] input input text like "href=hoge color=red pos="0.0 1.0""
 * @param[out] std::map A["href"]="hoge", A["pos"]="0.1 1.0"
 */
DFM2_INLINE void ParseAttributes
 (std::map<std::string, std::string>& mapAttr,
  const std::string& input);

/**
 * @brief separate tag and content in xml file
 * @details for xml text like "<h>hoge</h>", the tag is "h" and "/h" and content is "hoge".
 * tag is put at the odd index, and content is put at the even index.
 * The space at the beginning of the tag and content is removed
 */
DFM2_INLINE void XML_SeparateTagContent
 (std::vector<std::string>& aStr,
  const std::vector<char>& input);


/**
 * @brief check  if c is alphabet (e.g., A-Za-z) or not
 */
DFM2_INLINE bool isAlphabet(char c);
//DFM2_INLINE bool isNumber(char c);
DFM2_INLINE bool isAlphabetUpper(int iascii);
DFM2_INLINE bool isNumeric(int iascii);

}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/funcs.cpp"
#endif

#endif /* FUNCS_H */
