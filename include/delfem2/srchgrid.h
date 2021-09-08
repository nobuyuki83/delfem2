#ifndef DFM2_SRCHGRID_H
#define DFM2_SRCHGRID_H

#include <vector>
#include <algorithm> // for sort

namespace delfem2 {

template <class T>
void insertion_sort_stl(std::vector<T>& vec)
{
  for(auto it = vec.begin(); it != vec.end(); it++){
    auto const insertion_point = std::upper_bound(vec.begin(), it, *it);
    std::rotate(insertion_point, it, it+1);
  }
}

class SearchGrid {
public:
  class CGrid2Obj {
  public:
    bool operator<(const CGrid2Obj &rhs) const {
      return igrid < rhs.igrid;
    }

  public:
    unsigned int igrid;
    unsigned int iobj;
  };

public:
  std::vector<CGrid2Obj> aGrid2Obj;
  std::vector<unsigned int> aGrid2Obj_ind;
  double h = 1;
  double bbmin[3] = {0, 0, 0};
  double bbmax[3] = {1, 1, 1};
  unsigned int nx = 1, ny = 1, nz = 1;
public:
  void Initialize(
      const double bbmin_[3],
      const double bbmax_[3],
      const double h_,
      const unsigned int nobj) {
    h = h_;
    bbmin[0] = bbmin_[0];
    bbmin[1] = bbmin_[1];
    bbmin[2] = bbmin_[2];
    bbmax[0] = bbmax_[0];
    bbmax[1] = bbmax_[1];
    bbmax[2] = bbmax_[2];
    nx = static_cast<unsigned int>(ceil((bbmax[0] - bbmin[0]) / h));
    ny = static_cast<unsigned int>(ceil((bbmax[1] - bbmin[1]) / h));
    nz = static_cast<unsigned int>(ceil((bbmax[2] - bbmin[2]) / h));
    aGrid2Obj.resize(nobj);
  }

public:
  void GridIndex(
      int &ix0, 
	  int &iy0, 
	  int &iz0,
      const double p[3]) const
  {
    ix0 = static_cast<int>(floor((p[0] - bbmin[0]) / h));
    if (ix0 < 0) { ix0 = 0; } else if (ix0 >= int(nx)) { ix0 = int(nx - 1); }
    iy0 = static_cast<int>(floor((p[1] - bbmin[1]) / h));
    if (iy0 < 0) { iy0 = 0; } else if (iy0 >= int(ny)) { iy0 = int(ny - 1); }
    iz0 = static_cast<int>(floor((p[2] - bbmin[2]) / h));
    if (iz0 < 0) { iz0 = 0; } else if (iz0 >= int(nz)) { iz0 = int(nz - 1); }
  }

  unsigned int GetGridIndex(const double p[3]) const {
    int ix0, iy0, iz0;
    GridIndex(ix0, iy0, iz0, p);
    return iz0 * ny * nx + iy0 * nx + ix0;
  }

  /*
  void SetObject(unsigned int iobj, const double p[3]) {
    int ix0, iy0, iz0;
    GridIndex(ix0, iy0, iz0, p);
    aGrid2Obj[iobj].igrid = iz0 * ny * nx + iy0 * nx + ix0;
    aGrid2Obj[iobj].iobj = iobj;
  }
   */

  void PostProcess(bool is_initial) {
    if( is_initial ) {
      std::sort(aGrid2Obj.begin(), aGrid2Obj.end()); // quick sort
    }
    else{
      insertion_sort_stl(aGrid2Obj);
    }
    const unsigned int ng = nx * ny * nz;
    aGrid2Obj_ind.resize(ng + 1);
    aGrid2Obj_ind[0] = 0;
    {
      unsigned int i0 = 0;
      for (unsigned int ig = 0; ig < ng; ++ig) {
        while (aGrid2Obj[i0].igrid <= ig && i0 < aGrid2Obj.size()) { i0++; }
        aGrid2Obj_ind[ig + 1] = i0;
      }
    }
#ifndef NDEBUG
    for (unsigned int ig = 0; ig < ng; ++ig) {
      for (unsigned int ii0 = aGrid2Obj_ind[ig]; ii0 < aGrid2Obj_ind[ig + 1]; ++ii0) {
        assert(aGrid2Obj[ii0].igrid == ig);
      }
    }
#endif
  }

  void GetOneRingNeighbor(
      std::vector<unsigned int> &aIP,
      const double p[3]) const {
    int ix0, iy0, iz0;
    this->GridIndex(ix0, iy0, iz0, p);
    aIP.resize(0);
    for (int i = -1; i < 2; ++i) {
      for (int j = -1; j < 2; ++j) {
        for (int k = -1; k < 2; ++k) {
          const int ix1 = ix0 + i;
          if (ix1 < 0 || ix1 >= int(nx)) continue;
          const int iy1 = iy0 + j;
          if (iy1 < 0 || iy1 >= int(ny)) continue;
          const int iz1 = iz0 + k;
          if (iz1 < 0 || iz1 >= int(nz)) continue;
          const unsigned int ig1 = iz1 * ny * nx + iy1 * nx + ix1;
          for (unsigned int ii0 = aGrid2Obj_ind[ig1]; ii0 < aGrid2Obj_ind[ig1 + 1]; ++ii0) {
            aIP.push_back(aGrid2Obj[ii0].iobj);
          }
        }
      }
    }
  }
};

}

#endif