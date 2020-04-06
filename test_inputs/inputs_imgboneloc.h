//
//  inputs_imgboneloc.h
//  examples_glfwold_static
//
//  Created by Nobuyuki Umetani on 2020-04-06.
//

#ifndef inputs_imgboneloc_h
#define inputs_imgboneloc_h

#include "delfem2/vec2.h"

namespace delfem2 {

void BoneLocs_SmplUglysweater
 (std::string& name_img_in_input_test,
  double& scale,
  std::vector< std::pair<double,CVec2d> >& aBoneLoc)
{
  name_img_in_input_test = "uglysweater.jpg";
  scale = 0.00166667;
  aBoneLoc.push_back( std::make_pair( 2, CVec2d(258,777)) ); // right hip
  aBoneLoc.push_back( std::make_pair(16, CVec2d(390,360)) ); // left shoulder
  aBoneLoc.push_back( std::make_pair(17, CVec2d(201,360)) ); // right shoulder
  aBoneLoc.push_back( std::make_pair(18, CVec2d(420,533)) ); // left elbow
  aBoneLoc.push_back( std::make_pair(20, CVec2d(468,413)) ); // left wrist
  aBoneLoc.push_back( std::make_pair(21, CVec2d(192,709)) ); // right wrist
}

}


#endif /* inputs_imgboneloc_h */
