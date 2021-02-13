
#ifndef PMS_NORMDIST_H
#define PMS_NORMDIST_H

class CamTransform{
public:
  float K0inv[9];
  float R01[9];
  float t01[3];
  float K1[9];
};

#endif