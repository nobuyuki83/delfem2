
#ifndef DFM2_OPENGL_OLD_ADF3_H
#define DFM2_OPENGL_OLD_ADF3_H

#include "delfem2/dfm2_inline.h"

namespace delfem2::opengl{

void Draw_Wire(const delfem2::AdaptiveDistanceField3::CNode &n) {
  ::glLineWidth(1);
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_LINES);
  double hw_ = n.hw_;
  double cent_[3] = {n.cent_[0], n.cent_[1], n.cent_[2]};
  ::glVertex3d(cent_[0] - hw_, cent_[1] - hw_, cent_[2] - hw_);
  ::glVertex3d(cent_[0] + hw_, cent_[1] - hw_, cent_[2] - hw_);

  ::glVertex3d(cent_[0] + hw_, cent_[1] - hw_, cent_[2] - hw_);
  ::glVertex3d(cent_[0] + hw_, cent_[1] + hw_, cent_[2] - hw_);

  ::glVertex3d(cent_[0] + hw_, cent_[1] + hw_, cent_[2] - hw_);
  ::glVertex3d(cent_[0] - hw_, cent_[1] + hw_, cent_[2] - hw_);

  ::glVertex3d(cent_[0] - hw_, cent_[1] + hw_, cent_[2] - hw_);
  ::glVertex3d(cent_[0] - hw_, cent_[1] - hw_, cent_[2] - hw_);
  ////
  ::glVertex3d(cent_[0] - hw_, cent_[1] - hw_, cent_[2] + hw_);
  ::glVertex3d(cent_[0] + hw_, cent_[1] - hw_, cent_[2] + hw_);

  ::glVertex3d(cent_[0] + hw_, cent_[1] - hw_, cent_[2] + hw_);
  ::glVertex3d(cent_[0] + hw_, cent_[1] + hw_, cent_[2] + hw_);

  ::glVertex3d(cent_[0] + hw_, cent_[1] + hw_, cent_[2] + hw_);
  ::glVertex3d(cent_[0] - hw_, cent_[1] + hw_, cent_[2] + hw_);

  ::glVertex3d(cent_[0] - hw_, cent_[1] + hw_, cent_[2] + hw_);
  ::glVertex3d(cent_[0] - hw_, cent_[1] - hw_, cent_[2] + hw_);
  ////
  ::glVertex3d(cent_[0] - hw_, cent_[1] - hw_, cent_[2] - hw_);
  ::glVertex3d(cent_[0] - hw_, cent_[1] - hw_, cent_[2] + hw_);

  ::glVertex3d(cent_[0] + hw_, cent_[1] - hw_, cent_[2] - hw_);
  ::glVertex3d(cent_[0] + hw_, cent_[1] - hw_, cent_[2] + hw_);

  ::glVertex3d(cent_[0] + hw_, cent_[1] + hw_, cent_[2] - hw_);
  ::glVertex3d(cent_[0] + hw_, cent_[1] + hw_, cent_[2] + hw_);

  ::glVertex3d(cent_[0] - hw_, cent_[1] + hw_, cent_[2] - hw_);
  ::glVertex3d(cent_[0] - hw_, cent_[1] + hw_, cent_[2] + hw_);

  ::glEnd();
}

void DrawThisAndChild_Wire(
    const delfem2::AdaptiveDistanceField3::CNode &n,
    const std::vector<delfem2::AdaptiveDistanceField3::CNode> &aNo) {
  //      std::cout << "ichild " << ichilds_[0] << " " << ichilds_[1] << " " << ichilds_[2] << " " << ichilds_[3] << std::endl;
  if (n.ichilds_[0] == -1) {
    Draw_Wire(n);
    return;
  }
  DrawThisAndChild_Wire(aNo[n.ichilds_[0]], aNo);
  DrawThisAndChild_Wire(aNo[n.ichilds_[1]], aNo);
  DrawThisAndChild_Wire(aNo[n.ichilds_[2]], aNo);
  DrawThisAndChild_Wire(aNo[n.ichilds_[3]], aNo);
  DrawThisAndChild_Wire(aNo[n.ichilds_[4]], aNo);
  DrawThisAndChild_Wire(aNo[n.ichilds_[5]], aNo);
  DrawThisAndChild_Wire(aNo[n.ichilds_[6]], aNo);
  DrawThisAndChild_Wire(aNo[n.ichilds_[7]], aNo);
}

}

//#ifndef DFM2_STATIC_LIBRARY
//  #include "delfem2/opengl/old/gizmo.cpp"
//#endif

#endif /* gizmo_glold_h */
