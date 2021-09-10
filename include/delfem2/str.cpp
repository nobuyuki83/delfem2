/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/str.h"

#include <sstream>
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>

// ----------------------------------------

// probably std::stroi is safer to use but it is only for C++11
DFM2_INLINE int delfem2::myStoi(const std::string& str){
  char* e;
  long d = std::strtol(str.c_str(),&e,0);
  return (int)d;
}

DFM2_INLINE double delfem2::myStod(const std::string& str){
  char* e;
  double d = std::strtod(str.c_str(),&e);
  return d;
}



// -----------------------------------------------------------

DFM2_INLINE bool delfem2::isAlphabet(char c){
  if ( c >= 65 && c <= 90  ){ return true; }
  if(  c >= 97 && c <= 122 ){ return true; }
  return false;
}

/*
DFM2_INLINE bool delfem2::isNumber(char c){
  return c >= 46 && c <= 57;
}
 */

DFM2_INLINE bool delfem2::isAlphabetUpper(int iascii){
  if( iascii >= 65 && iascii <= 90 ) return true;
  return false;
}

DFM2_INLINE bool delfem2::isNumeric(int iascii){
  if(iascii >= 48 && iascii <= 57 ) return true;
  return false;
}

DFM2_INLINE void delfem2::Split(
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

DFM2_INLINE std::vector<std::string> delfem2::Split(
    const std::string& str,
    char delimiter)
{
  std::vector<std::string> aToken;
  Split(aToken,str,delimiter);
  return aToken;
}

DFM2_INLINE std::vector<std::string> delfem2::Split(
	const std::string& str,
	const std::string& del)
{
  std::vector<std::string> aToken;
  unsigned int imark = 0;
  bool is_del0 = true;
  for(unsigned int i=0;i<str.size();++i){
    bool is_del1 = false;
    for(char j : del){
      if( str[i] == j ){ is_del1 = true; break; }
    }
    if( !is_del0 && is_del1 ){ // just got delimitner after text
      aToken.emplace_back(str.data()+imark,str.data()+i );
    }
    if( is_del0 && !is_del1 ){ // sequence of delimitner
      imark = i;
    }
    is_del0 = is_del1;
  }
  if( !is_del0 ){ // just got delimitner after text
    aToken.emplace_back(str.data()+imark,str.data()+str.size() );
  }
  return aToken;
}

// "(a,b),c,(d,e)" - > "(a,b)", "c", "(d,e)"
DFM2_INLINE std::vector<std::string> delfem2::Split_Parentheses(
	const std::string& str,
	char delimiter,
	const std::string& par)
{
  std::vector<std::string> aToken;
  if( par.size() != 2 ){ aToken.push_back(str); return aToken; }
  char cs = par[0];
  char ce = par[1];
  //
  unsigned int is=0;
  int ilevel = 0;
  for(unsigned int ie=0;ie<str.size();++ie){
    if( ie == str.size()-1 ){
      aToken.emplace_back(str.data()+is,str.data()+ie+1 );
    }
    if( str[ie] == cs ){ ilevel++; }
    if( str[ie] == ce ){ ilevel--; }
    if( str[ie] == delimiter && ilevel == 0 ){
      aToken.emplace_back(str.data()+is,str.data()+ie );
      is = ie+1;
    }
  }
  return aToken;
}

DFM2_INLINE std::vector<std::string> delfem2::Split_Quote(
    const std::string& str,
    char delimiter,
    char quote)
{
  std::vector<std::string> aToken;
  unsigned int is=0;
  bool is_in = false;
  for(unsigned int ie=0;ie<str.size();++ie){
    if( ie == str.size()-1 ){
      aToken.emplace_back(str.data()+is,str.data()+ie+1 );
    }
    if( str[ie] == quote ){ is_in = !is_in; }
    if( str[ie] == delimiter && !is_in ){
      if( str[is] != delimiter ) { // skip the sequence of the delimiter
        aToken.emplace_back(str.data() + is, str.data() + ie);
      }
      is = ie+1;
    }
  }
  return aToken;
}

DFM2_INLINE std::string delfem2::Replace(
    const std::string& str,
    const char cf,
    const char ct)
{
  const size_t n = str.size();
  //
  std::string ss(str);
  for(unsigned int i=0;i<n;++i){
    if( ss[i] != cf ){ continue; }
    ss[i] = ct;
  }
  return ss;
}

DFM2_INLINE std::string delfem2::Remove(
	const std::string& str,
	const std::string& del)
{
  const size_t n = str.size();
  const size_t ndel = del.size();
  //
  std::string ss;
  ss.reserve(n);
  for(unsigned int i=0;i<n;++i){
    bool is_del = false;
    for(unsigned int idel=0;idel<ndel;++idel){
      if( str[i] != del[idel] ){ continue; }
      is_del = true;
      break;
    }
    if( is_del ) continue;
//    int i0 = (int)str[i];
//    char ich = (char)i0;
//    std::cout << i0 << " " << ich << std::endl;
    ss.push_back(str[i]);
  }
  return ss;
}

DFM2_INLINE std::string delfem2::RemoveSpace(
    const std::string& str)
{
  return Remove(str," ");
}

DFM2_INLINE std::string delfem2::RemoveBeginning(
    const std::string& str,
    const std::string& del)
{
  const size_t n = str.size();
  const size_t ndel = del.size();
  //
  int istat = 0;
  for(unsigned int i=0;i<n;++i){
    bool is_del = false;
    for(unsigned int idel=0;idel<ndel;++idel){
      if( str[i] != del[idel] ){ continue; }
      is_del = true;
      break;
    }
    if( is_del ) continue;
    istat = i;
    break;
  }
  return std::string(str.begin()+istat,str.end());
}

DFM2_INLINE std::string Remove_Quote(
    const std::string& str,
    char quat)
{
  const size_t n = str.size();
  {
    int nq = 0;
    for(unsigned int i=0;i<n;++i){
      if( str[i] == quat ){ ++nq; }
    }
    if( nq < 2 ){ return str;}
  }
  unsigned int istat = 0;
  for(;istat<n;++istat){
    if( str[istat] == quat ){ break; }
  }
  int iend = (int)n-1;
  for(;iend>=0;--iend){
    if( str[iend] == quat ){ break; }
  }
  return std::string(str.begin()+istat+1,str.begin()+iend);
}

// GetEnclosed with "()" --> "(a,(b,c),d)" -> a,(b,c),d
DFM2_INLINE std::string delfem2::Get_Parentheses(
	const std::string& str,
	const std::string& par)
{
  if( par.size() != 2 ){ return std::string(); }
  char cs = par[0];
  char ce = par[1];
  const size_t n = str.size();
  //
  int iss = -1;
  for(unsigned int i=0;i<n;++i){
    if( str[i] == cs ){ iss = i; break; }
  }
  //
  if( iss == -1 ){ return std::string(); }
  for(int i=int(n-1);i>=iss;--i){
    if( str[i] == ce ){
      std::string ss(str.begin()+iss+1,str.begin()+i);
      return ss;
    }
  }
  return std::string();
}

// ----------------------------------------------------------

// Read somehting like this {"command":"meshing_polygon","aXY_vertex":"0,0,0,30,30,30,30,0"}
DFM2_INLINE std::map<std::string, std::string>
delfem2::ReadDictionary_Json(const std::string& strIn)
{
  const std::string& buff = RemoveSpace(strIn);
  const size_t n =  buff.length();
  const char* p = buff.data();
  std::map<std::string, std::string> res;
  assert(p[0]=='{');
  assert(p[1]=='"');
  int ipos_begin = 1;
  int ipos_middle = -1;
  int ipos_end = -1;
  for(unsigned int i=1;i<n;++i){
    if( p[i]==':' && p[i-1]=='"' ){
      ipos_middle = i;
    }
    if( (p[i]==',' && p[i-1]=='"') || (i==n-1) ){
      ipos_end = i;
      assert(buff[ipos_begin]=='"');
      assert(buff[ipos_middle+1]=='"');
      assert(buff[ipos_middle-1]=='"');
      assert(buff[ipos_end-1]=='"');
      std::string key(buff.begin()+ipos_begin+1,buff.begin()+ipos_middle-1);
      std::string val(buff.begin()+ipos_middle+2,buff.begin()+ipos_end-1);
      res.insert( std::make_pair(key, val) );
      ipos_begin = ipos_end + 1;
    }
  }
  return res;
}

DFM2_INLINE void delfem2::ReadVector_CSV
(std::vector<double>& aVal,
 const std::string& strIn)
{
  const size_t n = strIn.length();
  const char* p = strIn.data();
  unsigned int ipos0 = 0;
  double val;
  for(unsigned int i=0;i<n;++i){
    if( p[i] == ',' ){
      std::string sval(p+ipos0,p+i);
//      sscanf(sval.c_str(),"%lf",&val);
      val = std::stod(sval);
      aVal.push_back(val);
      ipos0 = i+1;
    }
  }
  if( ipos0 < n ){
    std::string sval(p+ipos0,p+n);
//    sscanf(sval.c_str(),"%lf",&val);
    val = std::stod(sval);
    aVal.push_back(val);
  }
}

// ----------------------

DFM2_INLINE std::string delfem2::getCmdOptionStr
 (char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
  {
    return std::string(*itr);
  }
  return "";
}

DFM2_INLINE int delfem2::getCmdOptionInt
 (char ** begin, char ** end, const std::string & option, int ndef)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
  {
    return myStoi(*itr);
  }
  return ndef;
}

DFM2_INLINE bool delfem2::cmdOptionExists
 (char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

DFM2_INLINE double delfem2::getCmdOptionDouble
 (char ** begin, char ** end, const std::string & option, int ddef)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end){ return myStod(*itr); }
  return ddef;
}

DFM2_INLINE void delfem2::XML_SeparateTagContent
(std::vector<std::string>& aStr,
 const std::vector<char>& input)
{
  std::vector<char> buffer = input;
  char* s = buffer.data();
  char* mark = s;
  int state = 1;
  while (*s) {
    if (*s == '<' && state == 1) {
      // Start of a tag
      *s++ = '\0';
      aStr.emplace_back(mark);
      mark = s;
      state = 0;
    }
    else if (*s == '>' && state == 0 ) {       // Start of a content or new tag.
      *s++ = '\0';
      aStr.emplace_back(mark);
      mark = s;
      state = 1;
    }
    else {
      s++;
    }
  }
  for(auto & is : aStr){
    is = RemoveBeginning(is, " ");
  }
}


DFM2_INLINE void delfem2::ParseAttributes
(std::map<std::string, std::string>& mapAttr,
 const std::string& input)
{
  std::vector<std::string> aS = Split_Quote(input, ' ', '\"' );
  /*
  for(int is=0;is<aS.size();++is){
    std::cout << is << " " << aS[is] << std::endl;
  }
   */
  for(const auto & is : aS){
    std::vector<std::string> aS1 = Split(is, '=');
    if( aS1.size() != 2 ) continue;
    std::string s1 = Remove_Quote(aS1[1], '\"');
    mapAttr.insert( std::make_pair(aS1[0],s1) );
  }
}


DFM2_INLINE std::string delfem2::Str_SVGPolygon
(const std::vector<double>& aXY,
 double scale)
{
  /*
  double min_x,max_x, min_y,max_y;
  min_x = max_x = aXY[0];
  min_y = max_y = aXY[1];
  for(unsigned int ixy=0;ixy<aXY.size()/2;++ixy){
    if( aXY[ixy*2+0] < min_x ){ min_x = aXY[ixy*2+0]; }
    if( aXY[ixy*2+0] > max_x ){ max_x = aXY[ixy*2+0]; }
    if( aXY[ixy*2+1] < min_y ){ min_y = aXY[ixy*2+1]; }
    if( aXY[ixy*2+1] > max_y ){ max_y = aXY[ixy*2+1]; }
  }
  int w0 = (max_x-min_x)*scale;
  int h0 = (max_y-min_y)*scale;
  std::ostringstream oss;
  oss << "<?xml version=\"1.0\"?>\n";
  oss << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"";
  oss << w0;
  oss << "\" height=\"";
  oss << h0;
  oss << "\" viewBox=\"0 0 " << w0 << " " << h0;
  oss <<  "\">\n  <polygon points=\"";
  for(unsigned int ixy=0;ixy<aXY.size()/2;++ixy){
      double x0 = (aXY[ixy*2+0]-min_x)*scale;
      double y0 = (max_y-aXY[ixy*2+1])*scale;
      oss << x0 << "," << y0 << " ";
    }
  oss << "\" fill=\"blue\"></polygon>\n</svg>";
  return oss.str();
   */
  std::ostringstream oss;
  oss << "<?xml version=\"1.0\"?>\n";
  oss << "<svg xmlns=\"http://www.w3.org/2000/svg\">\n";
  oss << "<polygon points=\"";
  for(size_t ixy=0;ixy<aXY.size()/2;++ixy){
    const double x0 = +aXY[ixy*2+0]*scale;
    const double y0 = -aXY[ixy*2+1]*scale;
    oss << x0 << "," << y0 << " ";
  }
  oss << "\" fill=\"blue\"></polygon>\n</svg>";
  return oss.str();
}
