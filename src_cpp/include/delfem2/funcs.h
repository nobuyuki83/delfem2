/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef FUNCS_H
#define FUNCS_H

#include <fstream>
#include <map>
#include <vector>
#include <string>

// --------------------
// command related

std::string getCmdOptionStr(char ** begin, char ** end, const std::string & option);
int getCmdOptionInt(char ** begin, char ** end, const std::string & option, int ndef);
double getCmdOptionDouble(char ** begin, char ** end, const std::string & option, int ddef);
bool cmdOptionExists(char** begin, char** end, const std::string& option);

// -----------------
// Path related

bool isFileExists(const std::string& fpath);
std::string pathRemoveExtension(const std::string& fpath);
std::string pathGetExtension(const std::string& fpath);
std::string getPathDir(const std::string& fpath);

// ---------------
// File related

bool ReadParam(std::vector<float>& aPara,
               const std::string& fname);

template <typename T>
bool WriteParam
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

std::string LoadFile(const std::string& fname);

// ---------------------
// string handling

std::vector<std::string> Split(const std::string& str,
                               char delimiter);
void Split(std::vector<std::string>& aToken,
           const std::string& str,
           char delimiter);
std::vector<std::string> Split(const std::string& str,
                               const std::string& del);
std::vector<std::string> Split_Parentheses(const std::string& str,
                                           char delimiter,
                                           const std::string& par);
std::vector<std::string> Split_Quote(const std::string& str,
                                     char delimiter,
                                     char quote);
std::string Get_Parentheses(const std::string& str,
                        const std::string& par);
std::string Replace(const std::string& str,
                    const char cf, const char ct);
std::string Remove(const std::string& str,
                   const std::string& del);
std::string RemoveSpace(const std::string& str);
std::string RemoveBeginning(const std::string& str,
                            const std::string& del);

std::map<std::string, std::string> ReadDictionary(const std::string& path);
std::map<std::string, std::string> ReadDictionary_Json(const std::string& strIn);

void ReadVector_CSV(std::vector<double>& str, const std::string& strIn);

// ------------------------------
// Python realted funcs

std::map<std::string, std::string> ReadDictionary_Python(const std::string& buff);

bool LoadNumpy_2DimF(int& ndim0, int& ndim1, std::vector<float>& aData,
                     const std::string& path);
bool LoadNumpy_2DimD(int& ndim0, int& ndim1, std::vector<double>& aData,
                     const std::string& path);

bool LoadNumpy_1DimF(int& ndim0, std::vector<float>& aData,
                     const std::string& path);


// ---------------------
// SVG related function

std::string Str_SVGPolygon(const std::vector<double>& aXY,
                           double scale);


// ---------------------------------------------------

bool GetFileContents(std::vector<char>& aC,
                     const std::string& fpath);

/**
 * @brief make associative array
 * @param[in] input input text like "href=hoge color=red pos="0.0 1.0""
 * @param[out] std::map A["href"]="hoge", A["pos"]="0.1 1.0"
 */
void ParseAttributes(std::map<std::string, std::string>& mapAttr,
                     const std::string& input);

/**
 * @brief separate tag and content in xml file
 * @details for xml text like "<h>hoge</h>", the tag is "h" and "/h" and content is "hoge".
 * tag is put at the odd index, and content is put at the even index.
 * The space at the beginning of the tag and content is removed
 */
void XML_SeparateTagContent(std::vector<std::string>& aStr,
                            const std::vector<char>& input);


/**
 * @brief check  if c is alphabet (e.g., A-Za-z) or not
 */
bool isAlphabet(char c);
bool isNumber(char c);



#endif /* FUNCS_H */
