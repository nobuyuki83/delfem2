
/**
 * Find the fundamental matrix F from set of pairs of corresponding point coordinates
 * aP0P1 = (x0,y0,x1,y1)
 * (x1, y1, 1) * F * (x0, y0, 1) == 0
 * @param F fundamental matrix (row major)
 * @param aP0P1
 */
void Fundamental(
    float F[9],
    const std::vector<float>& aP0P1)
{
  const unsigned int npc = aP0P1.size()/4;
  const auto fnpc = float(npc);
  assert( aP0P1.size() == npc*4 );
  float c0[2] = {0.f, 0.f}, c1[2] = {0.f, 0.f};
  for(unsigned int ipc=0;ipc<npc;++ipc) {
    c0[0] += aP0P1[ipc*4+0];
    c0[1] += aP0P1[ipc*4+1];
    c1[0] += aP0P1[ipc*4+2];
    c1[1] += aP0P1[ipc*4+3];
  }
  c0[0] /= fnpc;
  c0[1] /= fnpc;
  c1[0] /= fnpc;
  c1[1] /= fnpc;
  float s0[2] = {0.f, 0.f}, s1[2] = {0.f, 0.f}; // scaling factor
  for(unsigned int ipc=0;ipc<npc;++ipc) {
    s0[0] += fabsf( c0[0] - aP0P1[ipc*4+0] );
    s0[1] += fabsf( c0[1] - aP0P1[ipc*4+1] );
    s1[0] += fabsf( c1[0] - aP0P1[ipc*4+2] );
    s1[1] += fabsf( c1[1] - aP0P1[ipc*4+3] );
  }
  s0[0] = fnpc/s0[0];
  s0[1] = fnpc/s0[1];
  s1[0] = fnpc/s1[0];
  s1[1] = fnpc/s1[1];
  Eigen::MatrixXd A(9,9);
  A.setZero();
  for(unsigned int ipc=0;ipc<npc;++ipc) {
    double x0 = (aP0P1[ipc*4+0] - c0[0])*s0[0];
    double y0 = (aP0P1[ipc*4+1] - c0[1])*s0[1];
    double x1 = (aP0P1[ipc*4+2] - c1[0])*s1[0];
    double y1 = (aP0P1[ipc*4+3] - c1[1])*s1[1];
    const double v[9] = {
        x1 * x0, x1 * y0, x1,
        y1 * x0, y1 * y0, y1,
        x0, y0, 1.0};
    for(int i=0;i<9;++i){
      for(int j=0;j<9;++j){
        A(i,j) += v[i]*v[j];
      }
    }
  }
  Eigen::Matrix3d F0;
  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A);
    const Eigen::VectorXd& u = eigensolver.eigenvectors().col(0);
    F0(0,0) = u(0); F0(0,1) = u(1);  F0(0,2) = u(2);
    F0(1,0) = u(3); F0(1,1) = u(4);  F0(1,2) = u(5);
    F0(2,0) = u(6); F0(2,1) = u(7);  F0(2,2) = u(8);
  }
  {
    Eigen::Matrix3d T0;
    {
      T0.setZero();
      T0(0, 0) = s0[0];
      T0(0, 2) = -c0[0] * s0[0];
      T0(1, 1) = s0[1];
      T0(1, 2) = -c0[1] * s0[1];
      T0(2, 2) = 1.0;
    }
    Eigen::Matrix3d T1;
    {
      T1.setZero();
      T1(0, 0) = s1[0];
      T1(0, 2) = -c1[0] * s1[0];
      T1(1, 1) = s1[1];
      T1(1, 2) = -c1[1] * s1[1];
      T1(2, 2) = 1.0;
    }
    F0 = T1.transpose()*(F0*T0);
  }
  {
    Eigen::JacobiSVD< Eigen::Matrix3d > svd(F0, Eigen::ComputeFullU | Eigen::ComputeFullV );
    Eigen::Matrix3d D = svd.singularValues().asDiagonal();
    D(2,2) = 0.0;
    F0 = svd.matrixU()*D*(svd.matrixV().transpose());
  }
  F0 /= F0(2,2);
  F[0*3+0] = F0(0, 0);
  F[0*3+1] = F0(0, 1);
  F[0*3+2] = F0(0, 2);
  F[1*3+0] = F0(1, 0);
  F[1*3+1] = F0(1, 1);
  F[1*3+2] = F0(1, 2);
  F[2*3+0] = F0(2, 0);
  F[2*3+1] = F0(2, 1);
  F[2*3+2] = F0(2, 2);
}