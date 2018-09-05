#ifndef ADF_H
#define ADF_H

#include <vector>

// virtual input class
class CInputAdaptiveDistanceField3D
{
public:
  virtual double Projection(double px, double py, double pz) const = 0;
};

class CAdaptiveDistanceField3D
{
public:
  CAdaptiveDistanceField3D();
  ~CAdaptiveDistanceField3D();
  void SetUp(const CInputAdaptiveDistanceField3D& ct, double bb[6]);
  void Draw() const;
  void SetFaceColor(double r, double g, double b){ color_[0] = r; color_[1] = g; color_[2] = b; }
  virtual double Projection
  (double px, double py, double pz,
   double n[3]) const;
  virtual bool IntersectionPoint
  (double p[3],
   const double org[3], const double dir[3]) const { return false; }
  virtual void GetMesh(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ, double elen) const{}
  /////
  void BuildIsoSurface_MarchingCube();
  void BuildMarchingCubeEdge();
  void SetShowCage(bool is_show){ this->is_show_cage = is_show; }
private:
  class CNode
  {
  public:
    CNode();
    CNode(const CNode& no);
    void SetCornerDist(const CInputAdaptiveDistanceField3D& ct);
    void Draw_Wire() const;
    void DrawThisAndChild_Wire(const std::vector<CNode>& aNo) const ;
    void MakeChildTree(const CInputAdaptiveDistanceField3D& ct, std::vector<CNode>& aNo, double min_hw, double max_hw);
    double FindDistNormal
    (double px, double py, double pz,
     double n[3],
     const std::vector<CNode>& aNo) const;
    void GenerateIsoSurface
    (std::vector<double>& aTri,
     const std::vector<CNode>& aNo) const;
  public:
    double cent_[3];
    double hw_;
    int ichilds_[8];
    double dists_[8];
  };
private:
  std::vector<CNode> aNode;
  double dist_min, dist_max;
  unsigned int nIsoTri_;
  double* aIsoTri_;
  double* aIsoEdge_;
  bool is_show_cage;
  double color_[3];
};

#endif /* adaptive_distance_field_h */
