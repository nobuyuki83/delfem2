/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/color.h"

template<typename T0, typename T1>
DFM2_INLINE void delfem2::GetRGB_HSV(
    T0 &r, T0 &g, T0 &b,
    T1 h, T1 s, T1 v) {
  r = static_cast<T0>(v);
  g = static_cast<T0>(v);
  b = static_cast<T0>(v);
  if (s > 0.0f) {
    h *= 6.0f;
    const int i = (int) h;
    const T1 f = h - static_cast<T1>(i);
    switch (i) {
      default:
      case 0: {
        g *= static_cast<T0>(1 - s * (1 - f));
        b *= static_cast<T0>(1 - s);
        break;
      }
      case 1: {
        r *= static_cast<T0>(1 - s * f);
        b *= static_cast<T0>(1 - s);
        break;
      }
      case 2: {
        r *= static_cast<T0>(1 - s);
        b *= static_cast<T0>(1 - s * (1 - f));
        break;
      }
      case 3: {
        r *= static_cast<T0>(1 - s);
        g *= static_cast<T0>(1 - s * f);
        break;
      }
      case 4: {
        r *= static_cast<T0>(1 - s * (1 - f));
        g *= static_cast<T0>(1 - s);
        break;
      }
      case 5: {
        g *= static_cast<T0>(1 - s);
        b *= static_cast<T0>(1 - s * f);
        break;
      }
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::GetRGB_HSV(
    float &r, float &g, float &b,
    float h, float s, float v);
template void delfem2::GetRGB_HSV(
    double &r, double &g, double &b,
    double h, double s, double v);
template void delfem2::GetRGB_HSV(
    float &r, float &g, float &b,
    double h, double s, double v);
#endif

// ------------------------------------------------------------

DFM2_INLINE void delfem2::interpolateColor(
    CColor &Cout,
    float r,
    const CColor &C0,
    const CColor &C1) {
  Cout.r = (1 - r) * C0.r + r * C1.r;
  Cout.g = (1 - r) * C0.g + r * C1.g;
  Cout.b = (1 - r) * C0.b + r * C1.b;
  Cout.a = (1 - r) * C0.a + r * C1.a;
}

DFM2_INLINE void delfem2::heatmap(
    double input,
    double *color) {
  if (0 <= input && input <= 0.25) {
    color[0] = 0.0;
    color[1] = input * 4.0;
    color[2] = 1.0;
  } else if (0.25 < input && input <= 0.5) {
    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 2.0 - input * 4.0;
  } else if (0.5 < input && input <= 0.75) {
    color[0] = input * 4 - 2.0;
    color[1] = 1.0;
    color[2] = 0.0;
  } else if (0.75 < input && input <= 1.0) {
    color[0] = 1.0;
    color[1] = 4.0 - input * 4.0;
    color[2] = 0.0;
  } else if (1.0 < input) {
    color[0] = 1.0;
    color[1] = 0.0;
    color[2] = 0.0;
  } else {
    color[0] = 0.0;
    color[1] = 0.0;
    color[2] = 1.0;
  }
}

DFM2_INLINE delfem2::CColor delfem2::getColor(
    double input,
    const std::vector<std::pair<double, CColor> > &colorMap) {
  if (colorMap.empty()) return CColor::Black();
  if (input < colorMap[0].first) {
    return colorMap[0].second;
  }
  for (int ic = 0; ic < (int) colorMap.size() - 1; ++ic) {
    double val0 = colorMap[ic].first;
    double val1 = colorMap[ic + 1].first;
    if (val0 <= input && input <= val1) {
      auto rp = (float) ((input - val0) / (val1 - val0));
      CColor color;
      delfem2::interpolateColor(color, rp, colorMap[ic].second, colorMap[ic + 1].second);
      return color;
    }
  }
  return colorMap[colorMap.size() - 1].second;
}

DFM2_INLINE void delfem2::ColorMap_BlueGrayRed(
    std::vector<std::pair<double, CColor> > &colorMap,
    float min,
    float max) {
  double diff = (max - min) * 0.25;
  colorMap.emplace_back(min + diff * 0, CColor(0.0f, 0.0f, 1.0f, 1.0f)); // blue
  colorMap.emplace_back(min + diff * 1, CColor(0.0f, 0.2f, 1.0f, 1.0f));
  colorMap.emplace_back(min + diff * 2, CColor(0.5f, 0.5f, 0.5f, 1.0f)); // gray
  colorMap.emplace_back(min + diff * 3, CColor(1.0f, 0.2f, 0.0f, 1.0f));
  colorMap.emplace_back(min + diff * 4, CColor(1.0f, 0.0f, 0.0f, 1.0f)); // red
}

DFM2_INLINE void delfem2::ColorMap_BlueCyanGreenYellowRed(
    std::vector<std::pair<double, CColor> > &colorMap,
    float min,
    float max,
    float alpha) {
  double diff = (max - min) * 0.25;
  colorMap.emplace_back(min + diff * 0, CColor(0.0f, 0.0f, 1.0f, alpha)); // blue
  colorMap.emplace_back(min + diff * 1, CColor(0.0f, 1.0f, 1.0f, alpha)); // cyan
  colorMap.emplace_back(min + diff * 2, CColor(0.0f, 1.0f, 0.0f, alpha)); // green
  colorMap.emplace_back(min + diff * 3, CColor(1.0f, 1.0f, 0.0f, alpha)); // yellow
  colorMap.emplace_back(min + diff * 4, CColor(1.0f, 0.0f, 0.0f, alpha)); // red
}

DFM2_INLINE void delfem2::ColorMap_RedYellowGreenCyanBlue(
    std::vector<std::pair<double, CColor> > &colorMap,
    float min,
    float max) {
  double diff = (max - min) * 0.25;
  colorMap.emplace_back(min + diff * 0, CColor(1.0f, 0.0f, 0.0f, 1.0f)); // red
  colorMap.emplace_back(min + diff * 1, CColor(1.0f, 1.0f, 0.0f, 1.0f)); // yellow
  colorMap.emplace_back(min + diff * 2, CColor(0.0f, 1.0f, 0.0f, 1.0f)); // green
  colorMap.emplace_back(min + diff * 3, CColor(0.0f, 1.0f, 1.0f, 1.0f)); // cyan
  colorMap.emplace_back(min + diff * 4, CColor(0.0f, 0.0f, 1.0f, 1.0f)); // blue
}

// ----------------------------------------------------

DFM2_INLINE void delfem2::Write_Ply_Tri2DMesh_HeightColor(
    const std::string &fname,
    const std::vector<int> &aTri1,
    const std::vector<double> &aXY1,
    const std::vector<double> &aVal,
    std::vector<std::pair<double, CColor> > &colorMap) {
  const size_t np = aXY1.size() / 2;
  const size_t ntri = aTri1.size() / 3;
  std::ofstream fout;
  fout.open(fname.c_str(), std::ios::out);
  fout << "ply" << std::endl;
  fout << "format ascii 1.0" << std::endl;
  fout << "element vertex " << np << std::endl;
  fout << "property float x" << std::endl;
  fout << "property float y" << std::endl;
  fout << "property float z" << std::endl;
  fout << "property uchar red" << std::endl;
  fout << "property uchar green" << std::endl;
  fout << "property uchar blue" << std::endl;
  fout << "element face " << ntri << std::endl;
  fout << "property list uchar int vertex_indices" << std::endl;
  fout << "end_header" << std::endl;
  for (unsigned int ip = 0; ip < np; ++ip) {
    double v = aVal[ip];
    CColor c = delfem2::getColor(v, colorMap);
    int cr, cg, cb;
    c.getRGBChar(cr, cg, cb);
    fout << aXY1[ip * 2 + 0] << " " << aXY1[ip * 2 + 1] << " " << v << " ";
    fout << cr << " " << cg << " " << cb << std::endl;
  }
  for (unsigned int itri = 0; itri < ntri; ++itri) {
    fout << "3 " << aTri1[itri * 3 + 0] << " " << aTri1[itri * 3 + 1] << " " << aTri1[itri * 3 + 2] << std::endl;
  }
}


