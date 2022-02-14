#ifndef DFM2_RIGIDBODY
#define DFM2_RIGIDBODY

#include <iostream>
#include <set>
#include <cstdio>

#include "delfem2/vec3.h"
#include "delfem2/mat3.h"


// class of rigid body
class CRigidBody
{
public:
  /*
  CRigidBody(){
    u = CVector3(0,0,0);
    R.SetIdentity();
  }
   */
  CRigidBody(double m, std::vector<double> xyz_cg){
    assert(xyz_cg.size()==3);
    this->m = m;
    cg.x = xyz_cg[0];
    cg.y = xyz_cg[1];
    cg.z = xyz_cg[2];
    u = delfem2::CVec3d(0,0,0);
    R.SetIdentity();
  }
  void addCP(std::vector<double> p){
    assert( aCP.size() == aCForce.size() );
    aCP.push_back({p[0],p[1],p[2]});
    aCForce.resize(aCForce.size()+1);
  }
public:
  delfem2::CVec3d cg; // the center of gravity.
  double m;     // mass of this plate
  std::vector<delfem2::CVec3d> aCP; // contact point
  std::vector< std::pair<delfem2::CVec3d,delfem2::CVec3d> > aExForce; // (position,magnitude) of external forces
  //
  delfem2::CVec3d u; // deformation
  delfem2::CMat3d R; // rotation
  std::vector<delfem2::CVec3d> aCForce;
};

class CJoint {
public:
  CJoint(int irb0, int irb1, std::vector<double> xyz_p){
    assert(xyz_p.size()==3);
    this->irb0 = irb0;
    this->irb1 = irb1;
    p.x = xyz_p[0];
    p.y = xyz_p[1];
    p.z = xyz_p[2];
  }
public:
  delfem2::CVec3d p; // position
  int irb0;    // id of rigid body
  int irb1;    // id of rigid body
  ////
  double jp0[4]; // joint position
  double jp1[4]; // joint position
  ////
  delfem2::CVec3d linear;
  delfem2::CVec3d torque;
};

      
class CRigidBodyAssembly_Static
{
public:
    
  CRigidBodyAssembly_Static();
  CRigidBodyAssembly_Static(const std::vector<CRigidBody>& aRB,
                            const std::vector<CJoint>& aJ);
  
  void AddRigidBody(const double centre_of_mass[3],
                    const double mass,
                    const std::vector<double>& contact_points);
  void AddJoint(const double position[3],
                const int body_index1,
                const int body_index2);
  void ClearProblem();
  
  
  // Setting A example Problem
  void SetExample();
  
  //solving methods
  void Solve(){
    Solve_InterPlane();
  }
  void ComputeForces();
  std::vector<double> MinMaxXYZ() const;

//  void PrintJointForce();
  
  //display
//  void Draw();
//  void DrawFloorGL();
  
private: // non-static private functions
  void Solve_InterPlane();
  void SolveOneIteration();
private:  // static functions
  static delfem2::CVec3d rand_vec(double s);
  static delfem2::CMat3d rand_rot();
public:
  //members
  std::vector<CRigidBody> aRigidBody; // array of rigid body
  std::vector<CJoint> aJoint; // array of joint
//  std::vector<CPlate> aPlate;
  
  delfem2::CVec3d n; // normal direction of floor (should be an unit vector)
  delfem2::CVec3d gravity; // gravity
  double cont_stiff; // contact_stiffness (insensitive)
  double trans_stiff; // joint_translation_stiffness (insensitive)
  double rot_stiff; // joint_rotation_stiffness (insensitive)
  
  int nitr;
  double damping_ratio;
  ////
  bool is_draw_force;
  bool is_draw_skeleton;
  bool is_draw_deformed;
  bool is_draw_section;
  bool is_draw_grid;
  bool is_draw_section_moment;
  double scale_force;
  double scale_torque;
  };




#endif
