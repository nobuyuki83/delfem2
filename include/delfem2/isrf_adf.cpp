/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/isrf_adf.h"

#include <iostream>
#include <math.h>
#include <vector>

namespace delfem2 {
namespace adf {

const int edgeTable[256] = {
    0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0};

const int triTable[256][16] =
    {
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
        {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
        {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
        {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
        {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
        {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
        {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
        {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
        {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
        {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
        {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
        {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
        {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
        {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
        {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
        {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
        {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
        {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
        {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
        {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
        {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
        {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
        {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
        {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
        {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
        {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
        {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
        {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
        {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
        {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
        {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
        {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
        {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
        {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
        {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
        {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
        {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
        {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
        {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
        {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
        {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
        {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
        {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
        {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
        {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
        {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
        {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
        {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
        {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
        {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
        {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
        {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
        {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
        {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
        {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
        {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
        {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
        {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
        {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
        {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
        {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
        {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
        {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
        {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
        {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
        {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
        {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
        {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
        {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
        {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
        {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
        {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
        {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
        {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
        {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
        {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
        {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
        {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
        {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
        {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
        {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
        {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
        {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
        {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
        {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
        {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
        {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
        {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
        {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
        {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
        {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
        {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
        {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
        {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
        {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
        {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
        {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
    };

const double phexflg[8][3] = {
    {-1, -1, -1},
    {+1, -1, -1},
    {+1, +1, -1},
    {-1, +1, -1},
    {-1, -1, +1},
    {+1, -1, +1},
    {+1, +1, +1},
    {-1, +1, +1},
};

DFM2_INLINE void VertexInterp
    (double pIntp[3],
     const double cent_[3], double hw_,
     unsigned int ind0, unsigned int ind1,
     double dist0, double dist1) {
  const double p0[3] = {
      cent_[0] + hw_ * phexflg[ind0][0],
      cent_[1] + hw_ * phexflg[ind0][1],
      cent_[2] + hw_ * phexflg[ind0][2]};
  const double p1[3] = {
      cent_[0] + hw_ * phexflg[ind1][0],
      cent_[1] + hw_ * phexflg[ind1][1],
      cent_[2] + hw_ * phexflg[ind1][2]};
  const double r0 = +dist1 / (dist1 - dist0);
  const double r1 = -dist0 / (dist1 - dist0);
  pIntp[0] = p0[0] * r0 + p1[0] * r1;
  pIntp[1] = p0[1] * r0 + p1[1] * r1;
  pIntp[2] = p0[2] * r0 + p1[2] * r1;
}

}
}

delfem2::CADF3::CADF3() {
  this->is_show_cage = false;
  //	this->is_show_cage = true;
  this->nIsoTri_ = 0;
  this->aIsoTri_ = 0;
  this->aIsoEdge_ = 0;
  color_[0] = 1.0;
  color_[1] = 1.0;
  color_[2] = 1.0;
}

delfem2::CADF3::~CADF3() {
  if (this->aIsoTri_ != 0) { delete[] this->aIsoTri_; }
  if (this->aIsoEdge_ != 0) { delete[] this->aIsoEdge_; }
}

void delfem2::CADF3::SetUp(
    const CInput_ADF3 &ct,
    double bb[6]) {
  aNode.reserve(1024 * 64);
  aNode.resize(1);
  CNode no;
  {
    no.cent_[0] = (bb[0] + bb[1]) * 0.5;
    no.cent_[1] = (bb[2] + bb[3]) * 0.5;
    no.cent_[2] = (bb[4] + bb[5]) * 0.5;
    no.hw_ = (bb[1] - bb[0]) > (bb[3] - bb[2]) ? (bb[1] - bb[0]) * 0.5 : (bb[3] - bb[2]) * 0.5;
    no.hw_ = no.hw_ > (bb[5] - bb[4]) * 0.5 ? no.hw_ : (bb[5] - bb[4]) * 0.5;
    no.hw_ *= 1.1234;
    no.SetCornerDist(ct);
    no.MakeChildTree(ct, aNode, no.hw_ * (0.99 / 64.0), no.hw_ * (1.01 / 4.0));
    //		no.MakeChildTree(ct,aNode,no.hw_*(0.99/128.0),no.hw_*(1.01/4.0));
    //		no.MakeChildTree(ct,aNode,no.hw_*(0.99/32.0),no.hw_*(1.01/4.0));
  }
  aNode[0] = no;
  std::cout << "ADF Oct-tree Node Size : " << aNode.size() << std::endl;
  ////
  dist_min = no.dists_[0];
  dist_max = dist_min;
  for (unsigned int ino = 0; ino < aNode.size(); ino++) {
    for (unsigned int i = 0; i < 8; i++) {
      const double dist = aNode[ino].dists_[i];
      dist_min = (dist < dist_min) ? dist : dist_min;
      dist_max = (dist > dist_max) ? dist : dist_max;
    }
  }
  // std::cout << "dist min max" << dist_min << " " << dist_max << std::endl;
  if (aIsoTri_ != 0) {
    delete aIsoTri_;
    aIsoTri_ = 0;
  }
  nIsoTri_ = 0;
}

// return penetration depth (inside is positive)
double delfem2::CADF3::Projection(
    double px, double py, double pz,
    double n[3]) const // normal outward
{
  const CNode &no = aNode[0];
  if (fabs(px - no.cent_[0]) > no.hw_ || fabs(py - no.cent_[1]) > no.hw_ || fabs(pz - no.cent_[2]) > no.hw_) {
    n[0] = no.cent_[0] - px;
    n[1] = no.cent_[1] - py;
    n[2] = no.cent_[2] - pz;
    const double invlen = 1.0 / sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (unsigned int i = 0; i < 3; i++) { n[0] *= invlen; }
    return -no.hw_;
  }
  return no.FindDistNormal(px, py, pz, n, aNode);
}

void delfem2::CADF3::BuildIsoSurface_MarchingCube() {
  std::vector<double> aTri;
  aTri.reserve(1024 * 32);
  aNode[0].GenerateIsoSurface(aTri, aNode);
  if (this->aIsoTri_ != 0) delete aIsoTri_;
  aIsoTri_ = new double[aTri.size()];
  for (unsigned int i = 0; i < aTri.size(); i++) { aIsoTri_[i] = aTri[i]; }
  nIsoTri_ = (unsigned int) aTri.size() / 9;
}

void delfem2::CADF3::BuildMarchingCubeEdge() {
  if (nIsoTri_ == 0) { this->BuildIsoSurface_MarchingCube(); }
  if (this->aIsoEdge_ != 0) { delete aIsoEdge_; }
  aIsoEdge_ = new double[nIsoTri_ * 18];
  for (unsigned int itri = 0; itri < nIsoTri_; itri++) {
    for (unsigned int i = 0; i < 3; i++) {
      aIsoEdge_[itri * 18 + 0 + i] = aIsoTri_[itri * 9 + 0 + i];
      aIsoEdge_[itri * 18 + 3 + i] = aIsoTri_[itri * 9 + 3 + i];

      aIsoEdge_[itri * 18 + 6 + i] = aIsoTri_[itri * 9 + 3 + i];
      aIsoEdge_[itri * 18 + 9 + i] = aIsoTri_[itri * 9 + 6 + i];

      aIsoEdge_[itri * 18 + 12 + i] = aIsoTri_[itri * 9 + 6 + i];
      aIsoEdge_[itri * 18 + 15 + i] = aIsoTri_[itri * 9 + 0 + i];
    }
  }
}

delfem2::CADF3::CNode::CNode() {
  cent_[0] = 0;
  cent_[1] = 0;
  cent_[2] = 0;
  hw_ = 0;
  ichilds_[0] = -1;
  ichilds_[1] = -1;
  ichilds_[2] = -1;
  ichilds_[3] = -1;
  ichilds_[4] = -1;
  ichilds_[5] = -1;
  ichilds_[6] = -1;
  ichilds_[7] = -1;
  dists_[0] = 0;
  dists_[1] = 0;
  dists_[2] = 0;
  dists_[3] = 0;
  dists_[4] = 0;
  dists_[5] = 0;
  dists_[6] = 0;
  dists_[7] = 0;
}

delfem2::CADF3::CNode::CNode(
    const CNode &no) {
  cent_[0] = no.cent_[0];
  cent_[1] = no.cent_[1];
  cent_[2] = no.cent_[2];
  hw_ = no.hw_;
  ////
  ichilds_[0] = no.ichilds_[0];
  ichilds_[1] = no.ichilds_[1];
  ichilds_[2] = no.ichilds_[2];
  ichilds_[3] = no.ichilds_[3];
  ichilds_[4] = no.ichilds_[4];
  ichilds_[5] = no.ichilds_[5];
  ichilds_[6] = no.ichilds_[6];
  ichilds_[7] = no.ichilds_[7];
  ////
  dists_[0] = no.dists_[0];
  dists_[1] = no.dists_[1];
  dists_[2] = no.dists_[2];
  dists_[3] = no.dists_[3];
  dists_[4] = no.dists_[4];
  dists_[5] = no.dists_[5];
  dists_[6] = no.dists_[6];
  dists_[7] = no.dists_[7];
}

void delfem2::CADF3::CNode::SetCornerDist(
    const CInput_ADF3 &ct) {
  dists_[0] = ct.sdf(cent_[0] - hw_, cent_[1] - hw_, cent_[2] - hw_);
  dists_[1] = ct.sdf(cent_[0] + hw_, cent_[1] - hw_, cent_[2] - hw_);
  dists_[2] = ct.sdf(cent_[0] + hw_, cent_[1] + hw_, cent_[2] - hw_);
  dists_[3] = ct.sdf(cent_[0] - hw_, cent_[1] + hw_, cent_[2] - hw_);
  dists_[4] = ct.sdf(cent_[0] - hw_, cent_[1] - hw_, cent_[2] + hw_);
  dists_[5] = ct.sdf(cent_[0] + hw_, cent_[1] - hw_, cent_[2] + hw_);
  dists_[6] = ct.sdf(cent_[0] + hw_, cent_[1] + hw_, cent_[2] + hw_);
  dists_[7] = ct.sdf(cent_[0] - hw_, cent_[1] + hw_, cent_[2] + hw_);
}

void delfem2::CADF3::CNode::MakeChildTree
    (const CInput_ADF3 &ct,
     std::vector<CNode> &aNo,
     double min_hw, double max_hw) {
  if (hw_ * 0.5 < min_hw) {
    ichilds_[0] = -1;
    return;
  }
  ////Edges
  const double va100 = ct.sdf(cent_[0], cent_[1] - hw_, cent_[2] - hw_);
  const double va210 = ct.sdf(cent_[0] + hw_, cent_[1], cent_[2] - hw_);
  const double va120 = ct.sdf(cent_[0], cent_[1] + hw_, cent_[2] - hw_);
  const double va010 = ct.sdf(cent_[0] - hw_, cent_[1], cent_[2] - hw_);

  const double va001 = ct.sdf(cent_[0] - hw_, cent_[1] - hw_, cent_[2]);
  const double va201 = ct.sdf(cent_[0] + hw_, cent_[1] - hw_, cent_[2]);
  const double va221 = ct.sdf(cent_[0] + hw_, cent_[1] + hw_, cent_[2]);
  const double va021 = ct.sdf(cent_[0] - hw_, cent_[1] + hw_, cent_[2]);

  const double va102 = ct.sdf(cent_[0], cent_[1] - hw_, cent_[2] + hw_);
  const double va212 = ct.sdf(cent_[0] + hw_, cent_[1], cent_[2] + hw_);
  const double va122 = ct.sdf(cent_[0], cent_[1] + hw_, cent_[2] + hw_);
  const double va012 = ct.sdf(cent_[0] - hw_, cent_[1], cent_[2] + hw_);

  ////Faces
  const double va101 = ct.sdf(cent_[0], cent_[1] - hw_, cent_[2]);
  const double va211 = ct.sdf(cent_[0] + hw_, cent_[1], cent_[2]);
  const double va121 = ct.sdf(cent_[0], cent_[1] + hw_, cent_[2]);
  const double va011 = ct.sdf(cent_[0] - hw_, cent_[1], cent_[2]);
  const double va110 = ct.sdf(cent_[0], cent_[1], cent_[2] - hw_);
  const double va112 = ct.sdf(cent_[0], cent_[1], cent_[2] + hw_);

  ////Center
  const double va111 = ct.sdf(cent_[0], cent_[1], cent_[2]);

  if (hw_ * 0.5 > max_hw) goto MAKE_CHILDS;

  double min_dist;
  {
    min_dist = fabs(va111);
    min_dist = (fabs(va100) < min_dist) ? fabs(va100) : min_dist;
    min_dist = (fabs(va210) < min_dist) ? fabs(va210) : min_dist;
    min_dist = (fabs(va120) < min_dist) ? fabs(va120) : min_dist;
    min_dist = (fabs(va010) < min_dist) ? fabs(va010) : min_dist;

    min_dist = (fabs(va001) < min_dist) ? fabs(va001) : min_dist;
    min_dist = (fabs(va201) < min_dist) ? fabs(va201) : min_dist;
    min_dist = (fabs(va221) < min_dist) ? fabs(va221) : min_dist;
    min_dist = (fabs(va021) < min_dist) ? fabs(va021) : min_dist;

    min_dist = (fabs(va101) < min_dist) ? fabs(va101) : min_dist;
    min_dist = (fabs(va211) < min_dist) ? fabs(va211) : min_dist;
    min_dist = (fabs(va121) < min_dist) ? fabs(va121) : min_dist;
    min_dist = (fabs(va011) < min_dist) ? fabs(va011) : min_dist;

    min_dist = (fabs(va102) < min_dist) ? fabs(va102) : min_dist;
    min_dist = (fabs(va212) < min_dist) ? fabs(va212) : min_dist;
    min_dist = (fabs(va122) < min_dist) ? fabs(va122) : min_dist;
    min_dist = (fabs(va012) < min_dist) ? fabs(va012) : min_dist;
    min_dist = (fabs(va110) < min_dist) ? fabs(va110) : min_dist;
    min_dist = (fabs(va112) < min_dist) ? fabs(va112) : min_dist;
  }

  if (min_dist > hw_ * 1.8) { // there is no mesh inside
    ichilds_[0] = -1;    // no-child
    return;
  }

  {
    if (min_dist < min_hw) goto MAKE_CHILDS;
  }

  {
    double t = min_hw * 0.8;
    if (fabs(va100 - (dists_[0] + dists_[1]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va210 - (dists_[1] + dists_[2]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va120 - (dists_[2] + dists_[3]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va010 - (dists_[3] + dists_[0]) * 0.5) > t) goto MAKE_CHILDS;

    if (fabs(va102 - (dists_[4] + dists_[5]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va212 - (dists_[5] + dists_[6]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va122 - (dists_[6] + dists_[7]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va012 - (dists_[7] + dists_[4]) * 0.5) > t) goto MAKE_CHILDS;

    if (fabs(va001 - (dists_[0] + dists_[4]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va201 - (dists_[1] + dists_[5]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va221 - (dists_[2] + dists_[6]) * 0.5) > t) goto MAKE_CHILDS;
    if (fabs(va021 - (dists_[3] + dists_[7]) * 0.5) > t) goto MAKE_CHILDS;

    if (fabs(va101 - (dists_[0] + dists_[1] + dists_[4] + dists_[5]) * 0.25) > t) goto MAKE_CHILDS;
    if (fabs(va211 - (dists_[1] + dists_[2] + dists_[5] + dists_[6]) * 0.25) > t) goto MAKE_CHILDS;
    if (fabs(va121 - (dists_[2] + dists_[3] + dists_[6] + dists_[7]) * 0.25) > t) goto MAKE_CHILDS;
    if (fabs(va011 - (dists_[3] + dists_[0] + dists_[7] + dists_[4]) * 0.25) > t) goto MAKE_CHILDS;
    if (fabs(va110 - (dists_[0] + dists_[1] + dists_[2] + dists_[3]) * 0.25) > t) goto MAKE_CHILDS;
    if (fabs(va112 - (dists_[4] + dists_[5] + dists_[6] + dists_[7]) * 0.25) > t) goto MAKE_CHILDS;

    if (fabs(
        va111 - (dists_[0] + dists_[1] + dists_[2] + dists_[3] + dists_[4] + dists_[5] + dists_[6] + dists_[7]) * 0.125)
        > t)
      goto MAKE_CHILDS;
  }

  ichilds_[0] = -1;    // no-child
  return;
  MAKE_CHILDS:
  const unsigned int nchild0 = static_cast<unsigned int>(aNo.size());
  aNo.resize(aNo.size() + 8);
  {    // left-bottom
    ichilds_[0] = nchild0;
    CNode no;
    no.cent_[0] = cent_[0] - hw_ * 0.5;
    no.cent_[1] = cent_[1] - hw_ * 0.5;
    no.cent_[2] = cent_[2] - hw_ * 0.5;
    no.hw_ = hw_ * 0.5;
    no.dists_[0] = dists_[0];
    no.dists_[1] = va100;
    no.dists_[2] = va110;
    no.dists_[3] = va010;
    no.dists_[4] = va001;
    no.dists_[5] = va101;
    no.dists_[6] = va111;
    no.dists_[7] = va011;
    no.MakeChildTree(ct, aNo, min_hw, max_hw);
    aNo[ichilds_[0]] = no;
  }
  {    // right-bottom
    ichilds_[1] = nchild0 + 1;
    CNode no;
    no.cent_[0] = cent_[0] + hw_ * 0.5;
    no.cent_[1] = cent_[1] - hw_ * 0.5;
    no.cent_[2] = cent_[2] - hw_ * 0.5;
    no.hw_ = hw_ * 0.5;
    no.dists_[0] = va100;
    no.dists_[1] = dists_[1];
    no.dists_[2] = va210;
    no.dists_[3] = va110;
    no.dists_[4] = va101;
    no.dists_[5] = va201;
    no.dists_[6] = va211;
    no.dists_[7] = va111;
    no.MakeChildTree(ct, aNo, min_hw, max_hw);
    aNo[ichilds_[1]] = no;
  }
  {    // right-top
    ichilds_[2] = nchild0 + 2;
    CNode no;
    no.cent_[0] = cent_[0] + hw_ * 0.5;
    no.cent_[1] = cent_[1] + hw_ * 0.5;
    no.cent_[2] = cent_[2] - hw_ * 0.5;
    no.hw_ = hw_ * 0.5;
    no.dists_[0] = va110;
    no.dists_[1] = va210;
    no.dists_[2] = dists_[2];
    no.dists_[3] = va120;
    no.dists_[4] = va111;
    no.dists_[5] = va211;
    no.dists_[6] = va221;
    no.dists_[7] = va121;
    no.MakeChildTree(ct, aNo, min_hw, max_hw);
    aNo[ichilds_[2]] = no;
  }
  {    // left-top
    ichilds_[3] = nchild0 + 3;
    CNode no;
    no.cent_[0] = cent_[0] - hw_ * 0.5;
    no.cent_[1] = cent_[1] + hw_ * 0.5;
    no.cent_[2] = cent_[2] - hw_ * 0.5;
    no.hw_ = hw_ * 0.5;
    no.dists_[0] = va010;
    no.dists_[1] = va110;
    no.dists_[2] = va120;
    no.dists_[3] = dists_[3];
    no.dists_[4] = va011;
    no.dists_[5] = va111;
    no.dists_[6] = va121;
    no.dists_[7] = va021;
    no.MakeChildTree(ct, aNo, min_hw, max_hw);
    aNo[ichilds_[3]] = no;
  }

  {    // left-bottom
    ichilds_[4] = nchild0 + 4;
    CNode no;
    no.cent_[0] = cent_[0] - hw_ * 0.5;
    no.cent_[1] = cent_[1] - hw_ * 0.5;
    no.cent_[2] = cent_[2] + hw_ * 0.5;
    no.hw_ = hw_ * 0.5;
    no.dists_[0] = va001;
    no.dists_[1] = va101;
    no.dists_[2] = va111;
    no.dists_[3] = va011;
    no.dists_[4] = dists_[4];
    no.dists_[5] = va102;
    no.dists_[6] = va112;
    no.dists_[7] = va012;
    no.MakeChildTree(ct, aNo, min_hw, max_hw);
    aNo[ichilds_[4]] = no;
  }
  {    // right-bottom
    ichilds_[5] = nchild0 + 5;
    CNode no;
    no.cent_[0] = cent_[0] + hw_ * 0.5;
    no.cent_[1] = cent_[1] - hw_ * 0.5;
    no.cent_[2] = cent_[2] + hw_ * 0.5;
    no.hw_ = hw_ * 0.5;
    no.dists_[0] = va101;
    no.dists_[1] = va201;
    no.dists_[2] = va211;
    no.dists_[3] = va111;
    no.dists_[4] = va102;
    no.dists_[5] = dists_[5];
    no.dists_[6] = va212;
    no.dists_[7] = va112;
    no.MakeChildTree(ct, aNo, min_hw, max_hw);
    aNo[ichilds_[5]] = no;
  }
  {    // right-top
    ichilds_[6] = nchild0 + 6;
    CNode no;
    no.cent_[0] = cent_[0] + hw_ * 0.5;
    no.cent_[1] = cent_[1] + hw_ * 0.5;
    no.cent_[2] = cent_[2] + hw_ * 0.5;
    no.hw_ = hw_ * 0.5;
    no.dists_[0] = va111;
    no.dists_[1] = va211;
    no.dists_[2] = va221;
    no.dists_[3] = va121;
    no.dists_[4] = va112;
    no.dists_[5] = va212;
    no.dists_[6] = dists_[6];
    no.dists_[7] = va122;
    no.MakeChildTree(ct, aNo, min_hw, max_hw);
    aNo[ichilds_[6]] = no;
  }
  {    // left-top
    ichilds_[7] = nchild0 + 7;
    CNode no;
    no.cent_[0] = cent_[0] - hw_ * 0.5;
    no.cent_[1] = cent_[1] + hw_ * 0.5;
    no.cent_[2] = cent_[2] + hw_ * 0.5;
    no.hw_ = hw_ * 0.5;
    no.dists_[0] = va011;
    no.dists_[1] = va111;
    no.dists_[2] = va121;
    no.dists_[3] = va021;
    no.dists_[4] = va012;
    no.dists_[5] = va112;
    no.dists_[6] = va122;
    no.dists_[7] = dists_[7];
    no.MakeChildTree(ct, aNo, min_hw, max_hw);
    aNo[ichilds_[7]] = no;
  }
  //			std::cout << "built" << ichilds_[0] << " " << ichilds_[1] << " " << ichilds_[2] << " " << ichilds_[3] << std::endl;
  return;
}

double delfem2::CADF3::CNode::FindDistNormal
    (double px, double py, double pz,
     double n[3],
     const std::vector<CNode> &aNo) const // normal outward
{
  if (fabs(px - cent_[0]) >= hw_
      || fabs(py - cent_[1]) >= hw_
      || fabs(pz - cent_[2]) >= hw_) {
    n[0] = px - cent_[0];
    n[1] = py - cent_[1];
    n[2] = pz - cent_[2];
    const double dist = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    const double inv_dist = 1.0 / dist;
    n[0] *= inv_dist;
    n[1] *= inv_dist;
    n[2] *= inv_dist;
    return -dist;
    //    return dist;
  }
  if (ichilds_[0] == -1) {
    const double rx = (px - cent_[0]) / hw_;
    const double ry = (py - cent_[1]) / hw_;
    const double rz = (pz - cent_[2]) / hw_;
    double dist =
        ((1 - rx) * (1 - ry) * (1 - rz) * dists_[0]
            + (1 + rx) * (1 - ry) * (1 - rz) * dists_[1]
            + (1 + rx) * (1 + ry) * (1 - rz) * dists_[2]
            + (1 - rx) * (1 + ry) * (1 - rz) * dists_[3]
            + (1 - rx) * (1 - ry) * (1 + rz) * dists_[4]
            + (1 + rx) * (1 - ry) * (1 + rz) * dists_[5]
            + (1 + rx) * (1 + ry) * (1 + rz) * dists_[6]
            + (1 - rx) * (1 + ry) * (1 + rz) * dists_[7]) * 0.125;
    ////
    n[0] =
        (-(1 - ry) * (1 - rz) * dists_[0]
            + (1 - ry) * (1 - rz) * dists_[1]
            + (1 + ry) * (1 - rz) * dists_[2]
            - (1 + ry) * (1 - rz) * dists_[3]
            - (1 - ry) * (1 + rz) * dists_[4]
            + (1 - ry) * (1 + rz) * dists_[5]
            + (1 + ry) * (1 + rz) * dists_[6]
            - (1 + ry) * (1 + rz) * dists_[7]);
    ////
    n[1] =
        (-(1 - rx) * (1 - rz) * dists_[0]
            - (1 + rx) * (1 - rz) * dists_[1]
            + (1 + rx) * (1 - rz) * dists_[2]
            + (1 - rx) * (1 - rz) * dists_[3]
            - (1 - rx) * (1 + rz) * dists_[4]
            - (1 + rx) * (1 + rz) * dists_[5]
            + (1 + rx) * (1 + rz) * dists_[6]
            + (1 - rx) * (1 + rz) * dists_[7]);
    ////
    n[2] =
        (-(1 - rx) * (1 - ry) * dists_[0]
            - (1 + rx) * (1 - ry) * dists_[1]
            - (1 + rx) * (1 + ry) * dists_[2]
            - (1 - rx) * (1 + ry) * dists_[3]
            + (1 - rx) * (1 - ry) * dists_[4]
            + (1 + rx) * (1 - ry) * dists_[5]
            + (1 + rx) * (1 + ry) * dists_[6]
            + (1 - rx) * (1 + ry) * dists_[7]);
    ////
    const double invlen = 1.0 / sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] *= -invlen;
    n[1] *= -invlen;
    n[2] *= -invlen;
    return dist;
  }
  if (px < cent_[0]) {
    if (py < cent_[1]) {
      if (pz < cent_[2]) { return aNo[ichilds_[0]].FindDistNormal(px, py, pz, n, aNo); }
      else { return aNo[ichilds_[4]].FindDistNormal(px, py, pz, n, aNo); }
    } else {
      if (pz < cent_[2]) { return aNo[ichilds_[3]].FindDistNormal(px, py, pz, n, aNo); }
      else { return aNo[ichilds_[7]].FindDistNormal(px, py, pz, n, aNo); }
    }
  } else {
    if (py < cent_[1]) {
      if (pz < cent_[2]) { return aNo[ichilds_[1]].FindDistNormal(px, py, pz, n, aNo); }
      else { return aNo[ichilds_[5]].FindDistNormal(px, py, pz, n, aNo); }
    } else {
      if (pz < cent_[2]) { return aNo[ichilds_[2]].FindDistNormal(px, py, pz, n, aNo); }
      else { return aNo[ichilds_[6]].FindDistNormal(px, py, pz, n, aNo); }
    }
  }
}

void delfem2::CADF3::CNode::GenerateIsoSurface
    (std::vector<double> &aTri,
     const std::vector<CNode> &aNo) const {
  if (this->ichilds_[0] != -1) {    // Evaluate Child Node
    for (unsigned int i = 0; i < 8; i++) { aNo[ichilds_[i]].GenerateIsoSurface(aTri, aNo); }
    return;
  }
  // Evaluate Dis Node
  int cubeindex = 0;
  {
    if (dists_[0] < 0) cubeindex |= 1;
    if (dists_[1] < 0) cubeindex |= 2;
    if (dists_[2] < 0) cubeindex |= 4;
    if (dists_[3] < 0) cubeindex |= 8;
    if (dists_[4] < 0) cubeindex |= 16;
    if (dists_[5] < 0) cubeindex |= 32;
    if (dists_[6] < 0) cubeindex |= 64;
    if (dists_[7] < 0) cubeindex |= 128;
  }
  //	std::cout << "Cube Index" << cubeindex << std::endl;
  double vertlist[12][3];
  if (adf::edgeTable[cubeindex] & 1) adf::VertexInterp(vertlist[0], cent_, hw_, 0, 1, dists_[0], dists_[1]);
  if (adf::edgeTable[cubeindex] & 2) adf::VertexInterp(vertlist[1], cent_, hw_, 1, 2, dists_[1], dists_[2]);
  if (adf::edgeTable[cubeindex] & 4) adf::VertexInterp(vertlist[2], cent_, hw_, 2, 3, dists_[2], dists_[3]);
  if (adf::edgeTable[cubeindex] & 8) adf::VertexInterp(vertlist[3], cent_, hw_, 3, 0, dists_[3], dists_[0]);
  if (adf::edgeTable[cubeindex] & 16) adf::VertexInterp(vertlist[4], cent_, hw_, 4, 5, dists_[4], dists_[5]);
  if (adf::edgeTable[cubeindex] & 32) adf::VertexInterp(vertlist[5], cent_, hw_, 5, 6, dists_[5], dists_[6]);
  if (adf::edgeTable[cubeindex] & 64) adf::VertexInterp(vertlist[6], cent_, hw_, 6, 7, dists_[6], dists_[7]);
  if (adf::edgeTable[cubeindex] & 128) adf::VertexInterp(vertlist[7], cent_, hw_, 7, 4, dists_[7], dists_[4]);
  if (adf::edgeTable[cubeindex] & 256) adf::VertexInterp(vertlist[8], cent_, hw_, 0, 4, dists_[0], dists_[4]);
  if (adf::edgeTable[cubeindex] & 512) adf::VertexInterp(vertlist[9], cent_, hw_, 1, 5, dists_[1], dists_[5]);
  if (adf::edgeTable[cubeindex] & 1024) adf::VertexInterp(vertlist[10], cent_, hw_, 2, 6, dists_[2], dists_[6]);
  if (adf::edgeTable[cubeindex] & 2048) adf::VertexInterp(vertlist[11], cent_, hw_, 3, 7, dists_[3], dists_[7]);

  for (unsigned int i = 0; adf::triTable[cubeindex][i] != -1; i += 3) {
    const unsigned int ind0 = static_cast<unsigned int>(aTri.size());
    aTri.resize(ind0 + 9);
    aTri[ind0 + 0] = vertlist[adf::triTable[cubeindex][i]][0];
    aTri[ind0 + 1] = vertlist[adf::triTable[cubeindex][i]][1];
    aTri[ind0 + 2] = vertlist[adf::triTable[cubeindex][i]][2];
    aTri[ind0 + 3] = vertlist[adf::triTable[cubeindex][i + 1]][0];
    aTri[ind0 + 4] = vertlist[adf::triTable[cubeindex][i + 1]][1];
    aTri[ind0 + 5] = vertlist[adf::triTable[cubeindex][i + 1]][2];
    aTri[ind0 + 6] = vertlist[adf::triTable[cubeindex][i + 2]][0];
    aTri[ind0 + 7] = vertlist[adf::triTable[cubeindex][i + 2]][1];
    aTri[ind0 + 8] = vertlist[adf::triTable[cubeindex][i + 2]][2];
  }
}
