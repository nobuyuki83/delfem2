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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

// "(a,b),c,(d,e)" - > "(a,b)", "c", "(d,e)"
std::vector<std::string> Split
(const std::string& str,
 char delimiter,
 const std::string& par)
{
  std::vector<std::string> aToken;
  if( par.size() != 2 ){ aToken.push_back(str); return aToken; }
  char cs = par[0];
  char ce = par[1];
  const int n = str.size();
  ////
  int is=0, ie;
  int ilevel = 0;
  for(ie=is+1;ie<n;++ie){
    if( ie == n-1 ){
      assert( ilevel == 0 );
      aToken.push_back( std::string(str.data()+is,str.data()+ie) );
    }
    else if( str[ie] == delimiter && ilevel == 0 ){
      aToken.push_back( std::string(str.data()+is,str.data()+ie) );
      is = ie+1;
    }
    else if( str[ie] == cs ){ ilevel++; }
    else if( str[ie] == ce ){ ilevel--; }
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

// GetEnclosed with "()" --> "(a,(b,c),d)" -> a,(b,c),d
std::string GetEnclosed
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

std::string Replace
(const std::string& str,
 const char cf,
 const char ct)
{
  const int n = str.size();
  ///
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
  std::string ss(str.begin()+istat,str.end());
  return ss;
}

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
  return 0;
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
  if( fin.is_open() ) return true;
  return false;
}

std::string pathRemoveExtension(const std::string& fpath)
{
  std::vector<std::string> aToken;
  Split(aToken, fpath, '.');
  std::string sRes;
  for(int it=0;it<(int)aToken.size()-1;++it){
    sRes = sRes + aToken[it];
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
std::map<std::string, std::string> ReadDictionary_Python(const std::string& strIn)
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
  char magic_string[6]; // 6 bytes (0x93NUMPY)
  unsigned char major_version; // 1 byte
  unsigned char minor_version; // 1 byte
  unsigned short header_len; // 2 bytes
};


bool LoadNumpy
(int& ndim0, int& ndim1,
 FILE* fp)
{
  NPY npy;
  fread(&npy, sizeof(npy), 1, fp);
  
  { // check magic string
    unsigned char sMagic[7] = {0x93,'N','U','M','P','Y'};
    if( memcmp(npy.magic_string, sMagic, 6 ) != 0 ){ return false; }
  }
  
  ndim0 = 0;
  ndim1 = 0;
  { // format
    char buff[256];
    fread(buff, 1, npy.header_len, fp);
    std::map<std::string, std::string> map0 = ReadDictionary_Python(std::string(buff));
    std::string str_shape = map0["'shape'"];
    str_shape = GetEnclosed(str_shape,"()");
    std::vector<std::string> aToken = Split(str_shape,',');
    if( aToken.size() != 2 ){ return false; }
    ndim0 = myStoi(aToken[0]);
    ndim1 = myStoi(aToken[1]);
  }
  return true;
}


bool LoadNumpy_2DimF
(int& ndim0, int& ndim1, std::vector<float>& aData,
 const std::string& path)
{
  FILE* fp = fopen(path.c_str(),"r");
  if( fp == NULL ) { return false; }
  LoadNumpy(ndim0, ndim1, fp);
  int size = ndim0*ndim1;
  aData.resize(size);
  fread(&aData[0], sizeof(float), size, fp);
  return true;
}

bool LoadNumpy_2DimD
(int& ndim0, int& ndim1, std::vector<double>& aData,
 const std::string& path)
{
  FILE* fp = fopen(path.c_str(),"r");
  if( fp == NULL ) { return false; }
  LoadNumpy(ndim0, ndim1, fp);
  int size = ndim0*ndim1;
  aData.resize(size);
  fread(&aData[0], sizeof(double), size, fp);
  return true;
}

bool LoadNumpy_1DimF
(int& ndim0, std::vector<float>& aData,
 const std::string& path)
{
  FILE* fp = fopen(path.c_str(),"r");
  if( fp == NULL ) { return false; }
  
  NPY npy;
  fread(&npy, sizeof(npy), 1, fp);
  
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
    fread(buff, 1, npy.header_len, fp);
    //    std::cout << buff << "###" << strlen(buff) << std::endl;
    std::map<std::string, std::string> map0 = ReadDictionary_Python(std::string(buff));
    //    for(std::map<std::string, std::string>::iterator itr=map0.begin();itr!=map0.end();++itr){
    //      std::cout << itr->first << " --> " << itr->second << std::endl;
    //    }
    std::string str_shape = map0["'shape'"];
    str_shape = GetEnclosed(str_shape,"()");
    std::vector<std::string> aToken = Split(str_shape,',');
    if( aToken.size() != 1 ){ return false; }
    ndim0 = myStoi(aToken[0]);
  }
  //    std::cout << ndim0 << " " << ndim1 << std::endl;
  
  //////
  double size = ndim0;
  aData.resize(size);
  //  double* aRes = (double*)malloc( sizeof(double)*size );
  fread(&aData[0], sizeof(float), size, fp);
  return true;
}
