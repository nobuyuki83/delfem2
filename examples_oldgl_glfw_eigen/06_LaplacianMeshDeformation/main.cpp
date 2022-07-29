#include <vector>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/msh_io_obj.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/opengl/old/funcs.h"

// ---------------------------------

void DrawTriangleMesh(
    const std::vector<unsigned int>& tri_vtx,
    const Eigen::MatrixXd& coord){
  ::glBegin(GL_TRIANGLES);
  for(unsigned int it=0;it<tri_vtx.size()/3;++it){
    const Eigen::Vector3d p0 = coord.row(tri_vtx[it*3+0]);
    const Eigen::Vector3d p1 = coord.row(tri_vtx[it*3+1]);
    const Eigen::Vector3d p2 = coord.row(tri_vtx[it*3+2]);
    ::glNormal3dv((p1-p0).cross(p2-p0).normalized().data());
    ::glVertex3dv(p0.data());
    ::glVertex3dv(p1.data());
    ::glVertex3dv(p2.data());
  }
  ::glEnd();
}

int main() {
  std::vector<unsigned int> tri_vtx;  // triangle index
  Eigen::MatrixXd vtx_xyz_ref;  // coordinate of vertices at reference configuration
  { // load bunny mesh
    std::vector<double> vtx_xyz;
    delfem2::Read_Obj3(
        vtx_xyz, tri_vtx,
        std::string(SOURCE_DIR) + "/../../test_inputs/bunny_1k.obj");
    delfem2::Normalize_Points3(vtx_xyz, 2.0);
    delfem2::Rotate_Points3(vtx_xyz, -M_PI * 0.5, 0., 0.);
    vtx_xyz_ref = Eigen::MatrixXd(vtx_xyz.size()/3, 3);
    for(int ip=0; ip<vtx_xyz_ref.rows(); ++ip){
      vtx_xyz_ref(ip, 0) = vtx_xyz[ip*3+0];
      vtx_xyz_ref(ip, 1) = vtx_xyz[ip*3+1];
      vtx_xyz_ref(ip, 2) = vtx_xyz[ip*3+2];
    }
  }
  const unsigned int num_vtx = vtx_xyz_ref.rows();  // number of vertices
  std::vector<unsigned int> psup, psup_ind; // jagged array for index of points surrounding points
  delfem2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      tri_vtx.data(), tri_vtx.size() / 3, 3, num_vtx);
  //
  std::vector<unsigned int> fix_base, fix_ear, fix_back; // index of fixed vertices
  Eigen::SparseMatrix<double> matrix_penalty(num_vtx, num_vtx); // sparse matrix for fix vertices constraint using penalty method
  { // make matrix for fixed constraint (diagonal is 1 if the corresponding vertex is fixed, otherwise zero)
    std::vector<Eigen::Triplet<double, int> > triplets;
    for(int ip=0; ip<vtx_xyz_ref.rows(); ++ip){
      if( vtx_xyz_ref(ip, 1) < -0.95 ){ // fixed vertices at base
        triplets.emplace_back(ip,ip,1);
        fix_base.push_back(ip);
      }
      if( vtx_xyz_ref(ip, 1) > +0.95 ){ // fixed vertices at ear
        triplets.emplace_back(ip,ip,1);
        fix_ear.push_back(ip);
      }
      if( vtx_xyz_ref(ip, 0) > +0.3 && vtx_xyz_ref(ip, 1) > +0.23 ){ // fixed vertices at back
        triplets.emplace_back(ip,ip,1);
        fix_back.push_back(ip);
      }
    }
    matrix_penalty.setFromTriplets(triplets.begin(), triplets.end()); // set sparse matrix's non-zero pattern and values
  }
  Eigen::SparseMatrix<double> matrix_laplacian(num_vtx, num_vtx); // graph laplacian matrix
  { // make laplacian matrix
    std::vector<Eigen::Triplet<double, int> > triplets;
    for(unsigned int ip=0;ip<psup_ind.size()-1;++ip){
      { // set diagonal value as valence of the vertex
        const auto val_dia = double(psup_ind[ip+1] - psup_ind[ip]);
        triplets.emplace_back(ip, ip, val_dia);
      }
      for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){ // set non diagonal vertex
        unsigned int jp0 = psup[ipsup]; // index of vertex connected to vertex "ip"
        triplets.emplace_back(ip, jp0, -1); // set non diagonal value as -1
      }
    }
    matrix_laplacian.setFromTriplets(triplets.begin(), triplets.end()); // set sparse matrix's non-zero pattern and values
  }
  const double penalty_coeff = 1.0e+5;
  const Eigen::SparseMatrix<double> sparse = matrix_laplacian.transpose() * matrix_laplacian + penalty_coeff * matrix_penalty;
  Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
  solver.analyzePattern(sparse);
  solver.factorize(sparse);
  const Eigen::MatrixXd Lv0 = matrix_laplacian * vtx_xyz_ref;

  // opengl starts here
  delfem2::glfw::CViewer3 viewer(2);
  //viewer.window_title = "task08";
  if (!glfwInit()) { exit(EXIT_FAILURE); }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  Eigen::MatrixXd vtx_xyz_def = vtx_xyz_ref; // deformed xyz coordinates of vertices
  while (!glfwWindowShouldClose(viewer.window)) {
    const double time = glfwGetTime();
    vtx_xyz_def = vtx_xyz_ref;
    for(unsigned int ip : fix_ear){ vtx_xyz_def(ip,0) += 0.5*sin(1*time+0); } // move points at ear
    for(unsigned int ip : fix_back){ vtx_xyz_def(ip,1) += 0.5*sin(2*time+0); } // move points at back
    Eigen::MatrixXd rhs0 = matrix_laplacian.transpose() * Lv0 + penalty_coeff * matrix_penalty * vtx_xyz_def;
    vtx_xyz_def = solver.solve(rhs0);
    // for each bone set relative rotation from parent bone
    // --------------------
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    delfem2::opengl::DrawAxis(1); // draw axes
    // draw fixed points
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    ::glColor3d(0,0,1); // draw fix points at base in blue
    for(unsigned int ip : fix_base){ ::glVertex3dv(Eigen::Vector3d(vtx_xyz_def.row(ip)).data()); }
    ::glColor3d(1,0,0);  // draw fix points at back in red
    for(unsigned int ip : fix_back){ ::glVertex3dv(Eigen::Vector3d(vtx_xyz_def.row(ip)).data()); }
    ::glColor3d(0,1,0);  // draw fix points at ear in green
    for(unsigned int ip : fix_ear){ ::glVertex3dv(Eigen::Vector3d(vtx_xyz_def.row(ip)).data()); }
    ::glEnd();
    // draw deformed mesh
    ::glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, std::array<float, 3>{1.f, 1.f, 1.f}.data());
    DrawTriangleMesh(tri_vtx,vtx_xyz_def);
    // draw end
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


