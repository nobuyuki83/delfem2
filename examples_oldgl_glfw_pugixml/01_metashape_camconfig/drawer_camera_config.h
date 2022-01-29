//
// Created by Nobuyuki Umetani on 2021-08-24.
//

#ifndef DRAWER_CAMERA_CONFIG_H
#define DRAWER_CAMERA_CONFIG_H

#include "delfem2/opengl/old/funcs.h"

void DrawCamera(
    const delfem2::MetashapeCamConfig::Camera& cam,
    const delfem2::MetashapeCamConfig& cam_config,
    float depth_scale)
{
  if( !cam.has_transform ){ return; }
//     delfem2::CMat4f mt = cam.transform.Inverse().transpose();
  delfem2::CMat4f mt = cam.transform.transpose();
  const unsigned int isnsr = cam.id_sensor;
  const auto& snsr = cam_config.aSensor[isnsr];
  const float z0 = depth_scale;
  const float x0 = (-0.5f*static_cast<float>(snsr.width)-snsr.cx)*z0/snsr.f;
  const float x1 = (+0.5f*static_cast<float>(snsr.width)-snsr.cx)*z0/snsr.f;
  const float y0 = (-0.5f*static_cast<float>(snsr.height)-snsr.cy)*z0/snsr.f;
  const float y1 = (+0.5f*static_cast<float>(snsr.height)-snsr.cy)*z0/snsr.f;
//      std::cout << "hoge" << snsr.width << " " << snsr.height << std::endl;
//      delfem2::CMat4f mt = cam.transform.transpose();
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ::glMultMatrixf(mt.mat);
  // ----
  // delfem2::opengl::DrawAxis(0.5);
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  ::glColor3d(1,0,0);
  ::glVertex3d(x0,y0,z0);  ::glVertex3d(x1,y0,z0);
  ::glColor3d(0,1,0);
  ::glVertex3d(x0,y0,z0);  ::glVertex3d(x0,y1,z0);
  ::glColor3d(0,0,0);
  ::glVertex3d(x1,y0,z0);  ::glVertex3d(x1,y1,z0);
  ::glVertex3d(x1,y1,z0);  ::glVertex3d(x0,y1,z0);
  ::glVertex3d(0,0,0); ::glVertex3d(x0,y0,z0);
  ::glVertex3d(0,0,0); ::glVertex3d(x1,y0,z0);
  ::glVertex3d(0,0,0); ::glVertex3d(x0,y1,z0);
  ::glVertex3d(0,0,0); ::glVertex3d(x1,y1,z0);
  ::glEnd();
  // ----
  ::glMatrixMode(GL_MODELVIEW);
  ::glPopMatrix();
}

void DrawRegion(
    const delfem2::MetashapeCamConfig& cam_config)
{
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  delfem2::CMat4f afm = delfem2::CMat4f::Mat3(cam_config.region_R.data());
  ::glTranslatef(cam_config.region_center.x,
                 cam_config.region_center.y,
                 cam_config.region_center.z);
  ::glMultMatrixf(afm.data());
  delfem2::opengl::DrawAABB3D_Edge(0,
                                   0,
                                   0,
                                   cam_config.region_size.x,
                                   cam_config.region_size.y,
                                   cam_config.region_size.z);
  ::glMatrixMode(GL_MODELVIEW);
  ::glPopMatrix();
}



#endif
