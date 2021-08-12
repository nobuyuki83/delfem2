/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_FILENPY_STR_H
#define DFM2_FILENPY_STR_H

#include <fstream>
#include <map>
#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

DFM2_INLINE std::map<std::string, std::string> ReadDictionary_Python(
    const std::string& buff);


template <typename REAL>
DFM2_INLINE bool LoadNumpy_2Dim(
    int& ndim0,
    int& ndim1,
    std::vector<REAL>& aData,
    const std::string& path);

DFM2_INLINE bool LoadNumpy_1DimF(
    int& ndim0,
    std::vector<float>& aData,
    const std::string& path);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/filenpy_str.cpp"
#endif

#endif /* FUNCS_H */
