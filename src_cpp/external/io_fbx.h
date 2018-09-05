#ifndef FBX_H
#define FBX_H

#include <vector>
#include <stdio.h>

#include "delfem2/rigmesh.h"

void Read_FBX(const std::string& path, CRigMsh& RigMsh);

#endif /* ioFBX_hpp */
