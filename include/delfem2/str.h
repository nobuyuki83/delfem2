/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file string processing (mainly for file parsing)
 */


#ifndef DFM2_STR_H
#define DFM2_STR_H

#include <map>
#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

DFM2_INLINE int myStoi(const std::string& str);
DFM2_INLINE double myStod(const std::string& str);

// --------------------
// command related

DFM2_INLINE std::string getCmdOptionStr(
    char ** begin,
    char ** end,
    const std::string & option);

DFM2_INLINE int getCmdOptionInt(
    char ** begin,
    char ** end,
    const std::string & option,
    int ndef);

DFM2_INLINE double getCmdOptionDouble(
    char ** begin,
    char ** end,
    const std::string & option,
    int ddef);

DFM2_INLINE bool cmdOptionExists(
    char** begin,
    char** end,
    const std::string& option);

// ---------------------
// string handling

DFM2_INLINE std::vector<std::string> Split(
    const std::string& str,
    char delimiter);

DFM2_INLINE void Split(
    std::vector<std::string>& aToken,
    const std::string& str,
    char delimiter);

DFM2_INLINE std::vector<std::string> Split(
    const std::string& str,
    const std::string& del);

DFM2_INLINE std::vector<std::string> Split_Parentheses(
    const std::string& str,
    char delimiter,
    const std::string& par);

// "'a,b',c,'d,e'" - > 'a,b' + 'c' + 'd,e' when delimitar is ',' and quote is '\''
DFM2_INLINE std::vector<std::string> Split_Quote(
    const std::string& str,
    char delimiter,
    char quote);

DFM2_INLINE std::string Get_Parentheses(
    const std::string& str,
    const std::string& par);

DFM2_INLINE std::string Replace(
    const std::string& str,
    const char cf, const char ct);

DFM2_INLINE std::string Remove(
    const std::string& str,
    const std::string& del);

DFM2_INLINE std::string RemoveSpace(
    const std::string& str);

DFM2_INLINE std::string RemoveBeginning(
    const std::string& str,
    const std::string& del);

DFM2_INLINE std::map<std::string, std::string> ReadDictionary_Json(
    const std::string& strIn);

DFM2_INLINE void ReadVector_CSV(
    std::vector<double>& str, const std::string& strIn);

// ---------------------
// SVG related function

DFM2_INLINE std::string Str_SVGPolygon(
    const std::vector<double>& aXY,
    double scale);


// ---------------------------------------------------
/**
 * @brief make associative array
 * @param[in] input input text like "href=hoge color=red pos="0.0 1.0""
 * @param[out] std::map A["href"]="hoge", A["pos"]="0.1 1.0"
 */
DFM2_INLINE void ParseAttributes(
    std::map<std::string, std::string>& mapAttr,
    const std::string& input);

/**
 * @brief separate tag and content in xml file
 * @details for xml text like "<h>hoge</h>", the tag is "h" and "/h" and content is "hoge".
 * tag is put at the odd index, and content is put at the even index.
 * The space at the beginning of the tag and content is removed
 */
DFM2_INLINE void XML_SeparateTagContent(
    std::vector<std::string>& aStr,
    const std::vector<char>& input);


/**
 * @brief check  if c is alphabet (e.g., A-Za-z) or not
 */
DFM2_INLINE bool isAlphabet(char c);
//DFM2_INLINE bool isNumber(char c);
DFM2_INLINE bool isAlphabetUpper(int iascii);
DFM2_INLINE bool isNumeric(int iascii);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/str.cpp"
#endif

#endif /* DFM2_STR_H */
