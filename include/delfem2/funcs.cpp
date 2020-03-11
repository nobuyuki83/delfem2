/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cassert>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include "delfem2/funcs.h"

namespace dfm2 = delfem2;

// ----------------------------------------

// probably std::stroi is safer to use but it is only for C++11
static int myStoi(const std::string& str){
  char* e;
  long d = std::strtol(str.c_str(),&e,0);
  return (int)d;
}

static double myStod(const std::string& str){
  char* e;
  double d = std::strtod(str.c_str(),&e);
  return d;
}



// -----------------------------------------------------------

bool isAlphabet(char c){
  if ( c >= 65 && c <= 90  ){ return true; }
  if(  c >= 97 && c <= 122 ){ return true; }
  return false;
}

bool isNumber(char c){
  return c >= 46 && c <= 57;
}

void Split
(std::vector<std::string>& aToken,
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

std::vector<std::string> Split
(const std::string& str,
 char delimiter)
{
  std::vector<std::string> aToken;
  Split(aToken,str,delimiter);
  return aToken;
}

std::vector<std::string> Split
(const std::string& str,
 const std::string& del)
{
  std::vector<std::string> aToken;
  int imark = 0;
  bool is_del0 = false;
  for(std::size_t i=0;i<str.size();++i){
    bool is_del1 = false;
    for(char j : del){
      if( str[i] == j ){ is_del1 = true; break; }
    }
    if( !is_del0 && is_del1 ){ // just got delimitner
      aToken.emplace_back(str.data()+imark,str.data()+i );
    }
    if( is_del0 && !is_del1 ){ // just got delimitner
      imark = i;
    }
    is_del0 = is_del1;
  }
  return aToken;
}

// "(a,b),c,(d,e)" - > "(a,b)", "c", "(d,e)"
std::vector<std::string> Split_Parentheses
(const std::string& str,
 char delimiter,
 const std::string& par)
{
  std::vector<std::string> aToken;
  if( par.size() != 2 ){ aToken.push_back(str); return aToken; }
  char cs = par[0];
  char ce = par[1];
  ////
  unsigned int is=0;
  int ilevel = 0;
  for(std::size_t ie=0;ie<str.size();++ie){
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

// "'a,b',c,'d,e'" - > 'a,b' + 'c' + 'd,e'
std::vector<std::string> Split_Quote
(const std::string& str,
 char delimiter,
 char quote)
{
  std::vector<std::string> aToken;
  unsigned int is=0;
  bool is_in = false;
  for(std::size_t ie=0;ie<str.size();++ie){
    if( ie == str.size()-1 ){
      aToken.emplace_back(str.data()+is,str.data()+ie+1 );
    }
    if( str[ie] == quote ){ is_in = !is_in; }
    if( str[ie] == delimiter && !is_in ){
      aToken.emplace_back(str.data()+is,str.data()+ie );
      is = ie+1;
    }
  }
  return aToken;
}




std::map<std::string, std::string> ReadDictionary(const std::string& fin_path)
{
  std::map<std::string, std::string> map0;
  ///
  std::ifstream fin;
  fin.open(fin_path.c_str());
  if( !fin.is_open() ){ return map0; }
  while(1){
    std::string str0, str1;
    fin >> str0 >> str1;
    if( fin.eof() ) break;
    //    std::cout << str0 << " " << str1 << std::endl;
    map0.insert( std::make_pair(str0,str1) );
  }
  return map0;
}

std::string Replace
(const std::string& str,
 const char cf,
 const char ct)
{
  const int n = str.size();
  //
  std::string ss(str);
  for(int i=0;i<n;++i){
    if( ss[i] != cf ){ continue; }
    ss[i] = ct;
  }
  return ss;
}

std::string Remove
(const std::string& str,
 const std::string& del)
{
  const int n = str.size();
  const int ndel = del.size();
  ///
  std::string ss;
  ss.reserve(n);
  for(int i=0;i<n;++i){
    bool is_del = false;
    for(int idel=0;idel<ndel;++idel){
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

std::string RemoveSpace
(const std::string& str)
{
  return Remove(str," ");
}

std::string RemoveBeginning
(const std::string& str,
 const std::string& del)
{
  const int n = str.size();
  const int ndel = del.size();
  ///
  int istat = 0;
  for(int i=0;i<n;++i){
    bool is_del = false;
    for(int idel=0;idel<ndel;++idel){
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

std::string Remove_Quote
(const std::string& str,
 char quat)
{
  const int n = str.size();
  {
    int nq = 0;
    for(int i=0;i<n;++i){
      if( str[i] == quat ){ ++nq; }
    }
    if( nq < 2 ){ return str;}
  }
  int istat = 0;
  for(;istat<n;++istat){
    if( str[istat] == quat ){ break; }
  }
  int iend = n-1;
  for(;iend>=0;--iend){
    if( str[iend] == quat ){ break; }
  }
  return std::string(str.begin()+istat+1,str.begin()+iend);
}

// GetEnclosed with "()" --> "(a,(b,c),d)" -> a,(b,c),d
std::string Get_Parentheses
(const std::string& str,
 const std::string& par)
{
  if( par.size() != 2 ){ return std::string(); }
  char cs = par[0];
  char ce = par[1];
  const int n = str.size();
  ////
  int iss = -1;
  for(int i=0;i<n;++i){
    if( str[i] == cs ){ iss = i; break; }
  }
  ////
  if( iss == -1 ){ return std::string(); }
  for(int i=n-1;i>=iss;--i){
    if( str[i] == ce ){
      std::string ss(str.begin()+iss+1,str.begin()+i);
      return ss;
    }
  }
  return std::string();
}

// ----------------------------------------------------------

// Read somehting like this {"command":"meshing_polygon","aXY_vertex":"0,0,0,30,30,30,30,0"}
std::map<std::string, std::string> ReadDictionary_Json(const std::string& strIn)
{
  const std::string& buff = RemoveSpace(strIn);
  const int n =  buff.length();
  const char* p = buff.data();
  std::map<std::string, std::string> res;
  assert(p[0]=='{');
  assert(p[1]=='"');
  int ipos_begin = 1;
  int ipos_middle = -1;
  int ipos_end = -1;
  for(int i=1;i<n;++i){
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

void ReadVector_CSV
(std::vector<double>& aVal,
 const std::string& strIn)
{
  const int n = strIn.length();
  const char* p = strIn.data();
  int ipos0 = 0;
  double val;
  for(int i=0;i<n;++i){
    if( p[i] == ',' ){
      std::string sval(p+ipos0,p+i);
      sscanf(sval.c_str(),"%lf",&val);
      aVal.push_back(val);
      ipos0 = i+1;
    }
  }
  if( ipos0 < n ){
    std::string sval(p+ipos0,p+n);
    sscanf(sval.c_str(),"%lf",&val);
    aVal.push_back(val);
  }
}

///////////////////////////////

std::string getCmdOptionStr(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
  {
    return std::string(*itr);
  }
  return "";
}

int getCmdOptionInt(char ** begin, char ** end, const std::string & option, int ndef)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
  {
    return myStoi(*itr);
  }
  return ndef;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

double getCmdOptionDouble(char ** begin, char ** end, const std::string & option, int ddef)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end){ return myStod(*itr); }
  return ddef;
}

//////////////////////////////////////////////////////////////////////////////
// file IO

std::string LoadFile
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

bool isFileExists(const std::string& fpath)
{
  std::ifstream fin;
  fin.open(fpath.c_str());
  return fin.is_open();
}

std::string pathRemoveExtension(const std::string& fpath)
{
  std::vector<std::string> aToken;
  Split(aToken, fpath, '.');
  std::string sRes;
  for(int it=0;it<(int)aToken.size()-1;++it){
    sRes += aToken[it];
  }
  return sRes;
}

std::string pathGetExtension(const std::string& fpath)
{
  std::vector<std::string> aToken;
  Split(aToken, fpath, '.');
  std::string sRes;
  if( !aToken.empty() ){
    sRes = aToken[aToken.size()-1];
  }
  return sRes;
}

std::string getPathDir(const std::string& fpath)
{
  const int iloc = fpath.find_last_of('/');
  std::string sres = std::string(fpath.begin(),fpath.begin()+iloc);
  return sres;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool ReadParam
(std::vector<float>& aPara,
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
    aPara.push_back(d);
  }
  return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Read somehting like this {'descr': '<f8', 'fortran_order': False, 'shape': (3682, 151), }
std::map<std::string, std::string>
dfm2::ReadDictionary_Python(const std::string& strIn)
{
  std::string buff = RemoveSpace(strIn);
  std::map<std::string, std::string> map0;
  /////
  const int n =  buff.length();
  const char* pbuff = buff.c_str();
  ////
  int iend = 1; // first letter is {
  int icolon = -1;
  bool is_parentheses = false;
  std::string skey;
  for(int i=0;i<n;++i){
    if( buff[i] == ':' ){
      char str[256]; strncpy(str, pbuff+iend,i-iend);
      str[i-iend] = '\0';
      //      std::cout << "key:" << str << "#" << std::endl;
      skey = std::string(str);
      //      if( buff[i+1] == ' ' ){ i+=1; }
      icolon = i;
    }
    if( buff[i] == ','){
      assert( icolon != -1 );
      if( is_parentheses ){ continue; }
      char str[256]; strncpy(str, pbuff+icolon+1,i-icolon-1);
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

bool isNaN(double x) { return x!=x; }

// 10bytes header
struct NPY
{
  char magic_string[6] = {'X','X','X','X','X','X'}; // 6 bytes (0x93NUMPY)
  unsigned char major_version = 0; // 1 byte
  unsigned char minor_version = 0; // 1 byte
  unsigned short header_len = 0; // 2 bytes
};


bool LoadNumpy
(int& ndim0, int& ndim1,
 FILE* fp)
{
  NPY npy;
  size_t n0 = fread(&npy, sizeof(npy), 1, fp);
  if( n0 != 1 ){ return false; }
  
  { // check magic string
    unsigned char sMagic[7] = {0x93,'N','U','M','P','Y'};
    if( memcmp(npy.magic_string, sMagic, 6 ) != 0 ){ return false; }
  }
  
  ndim0 = 0;
  ndim1 = 0;
  { // format
    char buff[256];
    size_t n1 = fread(buff, 1, npy.header_len, fp);
    if( n1 != npy.header_len ){ return false; }
    std::map<std::string, std::string> map0 = dfm2::ReadDictionary_Python(std::string(buff));
    std::string str_shape = map0["'shape'"];
    str_shape = Get_Parentheses(str_shape,"()");
    std::vector<std::string> aToken = Split(str_shape,',');
    if( aToken.size() != 2 ){ return false; }
    ndim0 = myStoi(aToken[0]);
    ndim1 = myStoi(aToken[1]);
  }
  return true;
}


bool dfm2::LoadNumpy_2DimF
(int& ndim0, int& ndim1, std::vector<float>& aData,
 const std::string& path)
{
  FILE* fp = fopen(path.c_str(),"rb");
  if( fp == nullptr ) { return false; }
  LoadNumpy(ndim0, ndim1, fp);
  int size = ndim0*ndim1;
  aData.resize(size);
  size_t n0 = fread(&aData[0], sizeof(float), size, fp);
  return (int) n0 == size;
}

bool dfm2::LoadNumpy_2DimD
(int& ndim0, int& ndim1, std::vector<double>& aData,
 const std::string& path)
{
  FILE* fp = fopen(path.c_str(),"rb");
  if( fp == nullptr ) { return false; }
  LoadNumpy(ndim0, ndim1, fp);
  int size = ndim0*ndim1;
  aData.resize(size);
  size_t n0 = fread(&aData[0], sizeof(double), size, fp);
  return (int) n0 == size;
}

bool dfm2::LoadNumpy_1DimF
(int& ndim0, std::vector<float>& aData,
 const std::string& path)
{
  FILE* fp = fopen(path.c_str(),"r");
  if( fp == nullptr ) { return false; }
  
  NPY npy;
  size_t n0 = fread(&npy, sizeof(npy), 1, fp);
  if( n0 != 1 ){ return false; }
  
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
    size_t n1 = fread(buff, 1, npy.header_len, fp);
    if( n1 != npy.header_len ){ return false; }
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
  double size = ndim0;
  aData.resize(size);
  //  double* aRes = (double*)malloc( sizeof(double)*size );
  size_t n2 = fread(&aData[0], sizeof(float), size, fp);
  return n2 == size;
}

// ----------------------

bool GetFileContents
(std::vector<char>& aC,
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


void XML_SeparateTagContent
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


void ParseAttributes
(std::map<std::string, std::string>& mapAttr,
 const std::string& input)
{
  std::vector<std::string> aS = Split_Quote(input, ' ', '\"' );
  for(const auto & is : aS){
    std::vector<std::string> aS1 = Split(is, '=');
    assert( aS1.size() == 2 );
    std::string s1 = Remove_Quote(aS1[1], '\"');
    mapAttr.insert( std::make_pair(aS1[0],s1) );
  }
}


std::string Str_SVGPolygon
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
    double x0 = +aXY[ixy*2+0]*scale;
    double y0 = -aXY[ixy*2+1]*scale;
    oss << x0 << "," << y0 << " ";
  }
  oss << "\" fill=\"blue\"></polygon>\n</svg>";
  return oss.str();
}
