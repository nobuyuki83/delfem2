#ifndef funcs_h
#define funcs_h

#include <map>

///////////////////////////////////////////////////
// command related

std::string getCmdOptionStr(char ** begin, char ** end, const std::string & option);
int getCmdOptionInt(char ** begin, char ** end, const std::string & option, int ndef);
double getCmdOptionDouble(char ** begin, char ** end, const std::string & option, int ddef);
bool cmdOptionExists(char** begin, char** end, const std::string& option);

///////////////////////////////////////////////////
// Path related

bool isFileExists(const std::string& fpath);
std::string pathRemoveExtension(const std::string& fpath);
std::string getPathDir(const std::string& fpath);

///////////////////////////////////////////////////
// File related

bool ReadParam(std::vector<float>& aPara,
               const std::string& fname);

bool WriteParam(const std::string& fname,
                const std::vector<float>& aPara);
std::string LoadFile(const std::string& fname);



//////////////////////////////////////////////////////////////////////////////
// string handling

std::vector<std::string> Split(const std::string& str,
                               char delimiter);
void Split(std::vector<std::string>& aToken,
           const std::string& str,
           char delimiter);
std::vector<std::string> Split(const std::string& str,
                               char delimiter,
                               const std::string& par);
std::string GetEnclosed(const std::string& str,
                        const std::string& par);
std::string Remove(const std::string& str,
                   const std::string& del);
std::string RemoveSpace(const std::string& str);

std::map<std::string, std::string> ReadDictionary(const std::string& path);
std::map<std::string, std::string> ReadDictionary_Json(const std::string& strIn);

void ReadVector_CSV(std::vector<double>& str, const std::string& strIn);

//////////////////////////////////////////////////////////////////////////////
// Python realted funcs

std::map<std::string, std::string> ReadDictionary_Python(const std::string& buff);

bool LoadNumpy_2DimF
(int& ndim0, int& ndim1, std::vector<float>& aData,
 const std::string& path);
bool LoadNumpy_2DimD
(int& ndim0, int& ndim1, std::vector<double>& aData,
 const std::string& path);

bool LoadNumpy_1DimF
(int& ndim0, std::vector<float>& aData,
 const std::string& path);


#endif /* fem_utility_h */
